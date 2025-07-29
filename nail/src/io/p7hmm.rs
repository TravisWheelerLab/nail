use std::{
    fs::File,
    io::{Read, Seek, SeekFrom},
    path::{Path, PathBuf},
    sync::Arc,
};

use crate::io::util::{ByteBufferExt, ReadSeekExt, ReadState, SeekableTake};

use libnail::structs::{
    profile::{ProfileBuilder, Transition},
    Profile,
};

use anyhow::{anyhow, bail, Context};
use indexmap::IndexMap;
use rayon::iter::{IntoParallelRefMutIterator, ParallelIterator};
use strum::{AsRefStr, EnumIter, EnumString, IntoEnumIterator};

pub fn profile_from_p7hmm_record_bytes(bytes: &[u8]) -> anyhow::Result<Profile> {
    enum ParseState {
        Header,
        Model,
    }

    #[derive(Debug)]
    enum ModelState {
        Comp,
        Mat,
        Ins,
        Trans,
    }

    let mut emission_buf = [0.0; Profile::MAX_ALPHABET_SIZE];
    let mut trans_buf = [0.0; Profile::NUM_STATE_TRANSITIONS - 1];
    let mut state = ParseState::Header;
    let mut model_state = ModelState::Comp;

    let mut prf = ProfileBuilder::default();
    let mut pos = 0;
    for line in bytes.split(|b| *b == b'\n').filter(|b| !b.is_empty()) {
        match state {
            ParseState::Header => {
                let flag = line
                    .first_word()?
                    .parse::<P7HeaderFlag>()
                    .unwrap_or_default();

                let value = std::str::from_utf8(&line[flag.as_ref().len()..])?.trim();

                match flag {
                    P7HeaderFlag::Name => {
                        prf.name(value.to_string());
                    }
                    P7HeaderFlag::Accession => {
                        prf.accession(value.to_string());
                    }
                    P7HeaderFlag::Length => {
                        prf.length(value.parse().with_context(|| {
                            format!("failed to parse HMM length: \"{}\"", value)
                        })?);
                    }
                    P7HeaderFlag::Alphabet => {
                        if value != "amino" {
                            bail!("non-amino HMM alphabet found: \"{}\"", value);
                        };
                    }
                    P7HeaderFlag::Stats => {
                        let stats: Vec<&str> = value.split_whitespace().collect();

                        if stats[1] == "FORWARD" {
                            prf.fwd_tau(stats[2].parse().with_context(|| {
                                format!("failed to parse HMM Forward tau: \"{}\"", stats[2])
                            })?);
                            prf.fwd_lambda(stats[3].parse().with_context(|| {
                                format!("failed to parse HMM Forward lambda: \"{}\"", stats[3])
                            })?);
                        }
                    }
                    P7HeaderFlag::Hmm => state = ParseState::Model,
                    _ => {
                        // ignore most flags for now
                    }
                }
            }
            ParseState::Model => {
                match model_state {
                    ModelState::Comp => {
                        if line.first_word()? == "COMPO" {
                            // note: we don't actually use the
                            //       composition in the p7HMM
                            //  why: HMMER doesn't
                            model_state = ModelState::Ins;
                        }
                    }
                    ModelState::Mat => {
                        let pos_pos = match line.first_non_whitespace_pos() {
                            Some(pos) => pos,
                            // line was empty
                            None => continue,
                        };
                        let pos_str = match line.word_from(pos_pos)? {
                            "//" => {
                                break;
                            }
                            s => s,
                        };

                        let new_pos = pos_str.parse::<usize>()?;

                        if new_pos != pos + 1 {
                            bail!("non-linear HMM positions: {pos}->{new_pos}");
                        } else {
                            pos = new_pos
                        }

                        parse_p7hmm_floats_into::<{ Profile::MAX_ALPHABET_SIZE }>(
                            std::str::from_utf8(&line[pos_pos + pos_str.len()..])?,
                            &mut emission_buf,
                        )?;

                        prf.mat_emission(pos, emission_buf);
                        model_state = ModelState::Ins;
                    }
                    ModelState::Ins => {
                        // note: we don't actually use the
                        //       insert emissions in the p7HMM
                        //  why: HMMER doesn't
                        model_state = ModelState::Trans;
                    }
                    ModelState::Trans => {
                        parse_p7hmm_floats_into::<{ Profile::NUM_STATE_TRANSITIONS - 1 }>(
                            std::str::from_utf8(line)?,
                            &mut trans_buf,
                        )?;

                        // map p7HMM transition order to Profile order
                        P7HmmTransition::iter().for_each(|t| {
                            trans_buf.swap(Transition::from(t) as usize, t as usize);
                        });

                        prf.transition(pos, trans_buf);

                        model_state = ModelState::Mat;
                    }
                }
            }
        }
    }

    prf.build()
}

#[derive(Default, PartialEq, EnumString, AsRefStr)]
#[strum(serialize_all = "UPPERCASE")]
enum P7HeaderFlag {
    #[default]
    #[strum(serialize = "")]
    None,
    #[strum(serialize = "HMMER3/f")]
    Format,
    Name,
    #[strum(serialize = "ACC")]
    Accession,
    #[strum(serialize = "DESC")]
    Description,
    #[strum(serialize = "LENG")]
    Length,
    #[strum(serialize = "MAXL")]
    MaxLength,
    #[strum(serialize = "ALPH")]
    Alphabet,
    #[strum(serialize = "RF")]
    Reference,
    #[strum(serialize = "MM")]
    Mask,
    #[strum(serialize = "CONS")]
    ConsensusResidue,
    #[strum(serialize = "CS")]
    ConsensusStructure,
    Map,
    Date,
    #[strum(serialize = "COM")]
    Command,
    Nseq,
    Effn,
    #[strum(serialize = "CKSUM")]
    Checksum,
    #[strum(serialize = "GA")]
    GatheringThresholds,
    #[strum(serialize = "TC")]
    TrustedCutoffs,
    #[strum(serialize = "NC")]
    NoiseCutoffs,
    Stats,
    Hmm,
}

#[repr(usize)]
#[derive(Clone, Copy, Debug, EnumIter)]
enum P7HmmTransition {
    MM = 0,
    MI = 1,
    MD = 2,
    IM = 3,
    II = 4,
    DM = 5,
    DD = 6,
}

impl From<P7HmmTransition> for Transition {
    fn from(p: P7HmmTransition) -> Self {
        use P7HmmTransition::*;
        match p {
            MM => Transition::MM,
            IM => Transition::IM,
            MI => Transition::MI,
            DM => Transition::DM,
            II => Transition::II,
            MD => Transition::MD,
            DD => Transition::DD,
        }
    }
}

fn negexp(x: f32) -> f32 {
    (-x).exp()
}

fn parse_p7hmm_floats_into<const N: usize>(floats: &str, out: &mut [f32; N]) -> anyhow::Result<()> {
    let mut tokens = floats.split_whitespace();

    out.iter_mut()
        .take(N)
        .try_for_each(|item| -> anyhow::Result<()> {
            *item = match tokens.next() {
                Some("*") => 0.0,
                Some(s) => s.parse().map(negexp)?,
                None => bail!("expected {N} floats, got fewer"),
            };

            Ok(())
        })
}

dyn_clone::clone_trait_object!(ProfileDatabase);
pub trait ProfileDatabase: dyn_clone::DynClone + Send + Sync {
    fn get(&mut self, name: &str) -> Option<Profile>;
    fn len(&self) -> usize;
    fn iter(&self) -> ProfileDatabaseIter;
}

pub struct ProfileDatabaseIter<'a> {
    pub(crate) inner: Box<dyn ProfileDatabase>,
    pub(crate) names_iter: Box<dyn DoubleEndedIterator<Item = &'a str> + 'a>,
}

impl<'a> ProfileDatabaseIter<'a> {
    pub fn names(self) -> Vec<&'a str> {
        self.names_iter.collect()
    }
}

impl<'a> Iterator for ProfileDatabaseIter<'a> {
    type Item = Profile;

    fn next(&mut self) -> Option<Self::Item> {
        self.names_iter.next().and_then(|n| self.inner.get(n))
    }

    // implementing size_hint to always return the
    // exact size of the iterator is REQUIRED for
    // the ExactSizerIterator trait to work; the
    // default impl of size_hint returns (0, None)
    fn size_hint(&self) -> (usize, Option<usize>) {
        let size = self.inner.len();
        (size, Some(size))
    }
}

// rayon expects the the iterators it
// uses to implement DoubleEndedIterator
impl<'a> DoubleEndedIterator for ProfileDatabaseIter<'a> {
    fn next_back(&mut self) -> Option<Self::Item> {
        self.names_iter.next_back().and_then(|n| self.inner.get(n))
    }
}

impl<'a> ExactSizeIterator for ProfileDatabaseIter<'a> {}

pub struct LexicalP7HmmIndex {
    pub(crate) offsets: IndexMap<String, P7HmmOffset>,
}

#[derive(Clone, Debug)]
pub struct P7HmmOffset {
    // points to the 'H' byte in the 'HMMER3/f ...' line
    start: usize,
    header_len_bytes: usize,
    model_len_bytes: usize,
    model_len: usize,
}

impl P7HmmOffset {
    pub fn new(start: usize) -> Self {
        Self {
            start,
            header_len_bytes: 0,
            model_len_bytes: 0,
            model_len: 0,
        }
    }
}

impl LexicalP7HmmIndex {
    fn from_path<P: AsRef<Path>>(path: P, n: usize) -> anyhow::Result<Self> {
        let mut file = File::open(path.as_ref())?;
        let sz = file.metadata().unwrap().len();
        let chunk_sz = sz / n as u64;

        let mut block_starts: Vec<u64> = (0..n).map(|i| i as u64).map(|i| i * chunk_sz).collect();
        // 2^16 gives us a ~65KiB buffer
        let mut buffer = vec![0; 2 << 16];

        for start in block_starts.iter_mut() {
            file.seek(SeekFrom::Start(*start)).expect("failed to seek");

            let mut offset = 0u64;
            'outer: loop {
                match file.read(&mut buffer) {
                    Ok(bytes_read) => {
                        for (pos, b) in buffer[0..bytes_read].iter().enumerate() {
                            // if we see an H, grab 8 characters
                            // and check if it's the format flag
                            if *b == b'H' {
                                if let Ok(P7HeaderFlag::Format) =
                                    buffer.as_slice().word_from(pos)?.parse::<P7HeaderFlag>()
                                {
                                    *start += offset;
                                    break 'outer;
                                }
                            }
                            offset += 1;
                        }
                    }
                    _ => bail!("failed to read from p7HMM buffer"),
                }
            }
        }

        let mut block_ends: Vec<u64> = block_starts.iter().skip(1).map(|b| b - 1).collect();
        block_ends.push(sz - 1);

        let mut files: Vec<_> = block_starts
            .into_iter()
            .zip(block_ends)
            .map(|(start, end)| {
                let file = File::open(path.as_ref())?;
                let limit = end - start + 1;
                Ok((SeekableTake::new(file, start, limit)?, start))
            })
            .collect::<anyhow::Result<_>>()?;

        let mut indexes = files
            .par_iter_mut()
            .map(|(file, start)| Self::new(file, Some(*start)))
            .collect::<anyhow::Result<Vec<_>>>()?;

        let mut index = indexes.remove(0);
        indexes.into_iter().for_each(|i| index.extend(i));

        index.offsets.shrink_to_fit();

        Ok(index)
    }

    fn new<R: Read + Seek>(mut data: R, start: Option<u64>) -> anyhow::Result<Self> {
        let mut buffer = [0; 8];

        match data.read(&mut buffer)? {
            8 => match std::str::from_utf8(&buffer)?.parse::<P7HeaderFlag>()? {
                P7HeaderFlag::Format => Ok(()),
                _ => Err(anyhow!("p7HMM buffer doesn't start with format flag")),
            },
            _ => Err(anyhow!("incomplete p7HMM buffer")),
        }
        .context("failed to validate p7HMM buffer")?;

        data.seek_relative(-8)?;

        #[derive(Debug)]
        enum ParseState {
            Format,
            Header,
            Model,
            Process,
        }

        let start = start.unwrap_or(0) as usize;
        let mut parse_state = ParseState::Header;
        let mut offsets = IndexMap::new();

        // 2<<20 == 2^21 | ~2MiB buffer
        let mut buffer = vec![0; 2 << 20];
        let mut name = String::new();
        let mut total_bytes_read = 1 + start;

        let mut offset = P7HmmOffset::new(start);

        let mut process_fn = |name: &str, offset| {
            if name.trim().is_empty() {
                bail!("failed to parse name from p7HMM buffer");
            }

            let mut name_clone = name.to_owned();
            name_clone.shrink_to_fit();
            offsets.insert(name_clone, offset);
            Ok(())
        };

        while let Ok(read_state) = data.read_with_state(&mut buffer) {
            let buf_slice = match read_state {
                ReadState::Reading(n) => {
                    let last_newline_pos =
                        buffer[..n]
                            .iter()
                            .rposition(|&b| b == b'\n')
                            .ok_or(anyhow!(
                                "failed to find any newlines in p7HMM buffer slice of length: {}",
                                n
                            ))?;

                    let n_bytes_back = (n - last_newline_pos) as i64;
                    data.seek_relative(-n_bytes_back)?;

                    // we don't want to include the last newline, or we'll
                    // likely end up triggering a flag parse on the void
                    &buffer[..last_newline_pos]
                }
                ReadState::Final(n) => &buffer[..n],
                ReadState::Done => {
                    match parse_state {
                        ParseState::Format => Ok(()),
                        s => Err(anyhow!("terminated p7HMM parsing in state:{:?} ", s)),
                    }?;
                    break;
                }
            };

            for (pos, &byte) in buf_slice.iter().enumerate() {
                let current_offset = total_bytes_read + pos;

                match parse_state {
                    ParseState::Format => {
                        if byte == b'\n' {
                            match buf_slice
                                .word_from(pos + 1)?
                                .parse::<P7HeaderFlag>()
                                .unwrap_or_default()
                            {
                                P7HeaderFlag::None => Ok(()),
                                P7HeaderFlag::Format => {
                                    parse_state = ParseState::Header;
                                    Ok(())
                                }
                                f => Err(anyhow!(
                                    "p7HMM record starts with incorrect flag: \"{}\"",
                                    f.as_ref()
                                )),
                            }?
                        }
                    }
                    ParseState::Header => {
                        if byte == b'\n' {
                            match buf_slice
                                .word_from(pos + 1)?
                                .parse::<P7HeaderFlag>()
                                .unwrap_or_default()
                            {
                                P7HeaderFlag::Name => {
                                    name.push_str(
                                        buf_slice.word_from(
                                            pos + P7HeaderFlag::Name.as_ref().len() + 3,
                                        )?,
                                    );
                                }
                                P7HeaderFlag::Length => {
                                    offset.model_len = buf_slice
                                        .word_from(pos + P7HeaderFlag::Length.as_ref().len() + 3)?
                                        .parse()
                                        .context("failed to parse length")?;
                                }
                                P7HeaderFlag::Hmm => {
                                    offset.header_len_bytes = current_offset - offset.start;
                                    parse_state = ParseState::Model;
                                }
                                _ => {}
                            }
                        }
                    }
                    ParseState::Model => {
                        if byte == b'/' {
                            if let "//" = buf_slice.word_from(pos)? {
                                offset.model_len_bytes =
                                    current_offset - (offset.start + offset.header_len_bytes) + 1;
                                parse_state = ParseState::Process
                            }
                        }
                    }
                    ParseState::Process => {
                        process_fn(&mut name, offset)?;

                        name.clear();
                        offset = P7HmmOffset::new(current_offset + 1);
                        parse_state = ParseState::Format;
                    }
                }
            }

            total_bytes_read += buf_slice.len();
        }

        offsets.shrink_to_fit();
        Ok(Self { offsets })
    }

    pub fn len(&self) -> usize {
        self.offsets.len()
    }

    fn offset_by_name(&self, name: &str) -> Option<P7HmmOffset> {
        self.offsets.get(name).cloned()
    }

    fn extend(&mut self, other: Self) {
        self.offsets.extend(other.offsets);
    }
}

pub struct P7Hmm {
    path: PathBuf,
    file: File,
    pub(crate) index: Arc<LexicalP7HmmIndex>,
    buffer: Vec<u8>,
}

impl Clone for P7Hmm {
    fn clone(&self) -> Self {
        let file = match File::open(&self.path) {
            Ok(file) => file,
            Err(err) => panic!(
                "failed to reopen fasta file on clone: {:?}\n error: {}",
                self.path, err
            ),
        };

        Self {
            file,
            path: self.path.clone(),
            index: self.index.clone(),
            buffer: vec![],
        }
    }
}

impl P7Hmm {
    pub fn from_path_par<P: AsRef<Path>>(path: P, n: usize) -> anyhow::Result<Self> {
        let index = Arc::new(LexicalP7HmmIndex::from_path(path.as_ref(), n)?);
        let file = File::open(path.as_ref())?;

        Ok(Self {
            file,
            index,
            buffer: Vec::new(),
            path: PathBuf::from(path.as_ref()),
        })
    }

    pub fn from_path<P: AsRef<Path>>(path: P) -> anyhow::Result<Self> {
        let mut file = File::open(path.as_ref())?;
        let index = Arc::new(LexicalP7HmmIndex::new(&mut file, None)?);

        Ok(Self {
            file,
            index,
            buffer: Vec::new(),
            path: PathBuf::from(path.as_ref()),
        })
    }

    pub fn len(&self) -> usize {
        self.index.len()
    }

    pub fn get(&mut self, name: &str) -> Option<Profile> {
        let offset = self.index.offset_by_name(name)?;

        self.buffer
            .resize(offset.model_len_bytes + offset.header_len_bytes, 0u8);

        self.file
            .seek(SeekFrom::Start(offset.start as u64))
            .expect("failed to seek in P7Hmm::get()");

        self.file
            .read_exact(&mut self.buffer)
            .expect("failed to read in P7Hmm::get()");

        Some(
            profile_from_p7hmm_record_bytes(&self.buffer).unwrap_or_else(|e| {
                panic!("failed to produce Profile in P7Hmm::get()\nError: {e}");
            }),
        )
    }

    pub fn names_iter(&self) -> impl Iterator<Item = &str> {
        self.index.offsets.keys().map(|k| k.as_str())
    }

    pub fn lengths_iter(&self) -> impl Iterator<Item = usize> + use<'_> {
        self.index.offsets.values().map(|v| v.model_len)
    }
}

impl ProfileDatabase for P7Hmm {
    fn get(&mut self, name: &str) -> Option<Profile> {
        self.get(name)
    }

    fn len(&self) -> usize {
        self.len()
    }

    fn iter(&self) -> ProfileDatabaseIter {
        ProfileDatabaseIter {
            inner: Box::new(self.clone()),
            names_iter: Box::new(self.index.offsets.keys().map(|s| s.as_str())),
        }
    }
}
