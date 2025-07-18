mod impl_rayon;

use std::{
    fs::File,
    io::{Read, Seek},
    path::{Path, PathBuf},
    sync::Arc,
};

use anyhow::{anyhow, bail, Context};
use indexmap::IndexMap;
use libnail::{
    alphabet::UTF8_TO_DIGITAL_AMINO,
    structs::{Profile, Sequence},
};
use rayon::iter::{IntoParallelRefMutIterator, ParallelIterator};
use strum_macros::{AsRefStr, EnumString};

fn profile_from_p7hmm_record_bytes(bytes: &[u8]) -> anyhow::Result<Profile> {
    enum ParseState {
        Header,
        Model,
    }

    enum ModelState {
        Comp,
        Mat,
        Ins,
        Trans,
    }

    #[derive(EnumString, AsRefStr)]
    #[strum(serialize_all = "UPPERCASE")]
    enum P7HeaderFlag {
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
        #[strum(serialize = "CSKSUM")]
        Ck,
        #[strum(serialize = "GA")]
        GatheringThreshold,
        #[strum(serialize = "TC")]
        TrustedCutoffs,
        #[strum(serialize = "NC")]
        NoiseCutoffs,
        Stats,
        Hmm,
    }

    fn unlog(x: f32) -> f32 {
        (-x).exp()
    }

    fn parse_floats_into<const N: usize>(floats: &str, out: &mut [f32; N]) -> anyhow::Result<()> {
        let mut tokens = floats.split_whitespace();

        out.iter_mut()
            .take(N)
            .try_for_each(|item| -> anyhow::Result<()> {
                *item = tokens
                    .next()
                    .ok_or_else(|| anyhow!("expected {N} floats, got fewer"))?
                    .parse()
                    .map(unlog)?;
                Ok(())
            })?;

        match tokens.next() {
            Some(_) => Err(anyhow!("expected {N} floats, got more")),
            None => Ok(()),
        }
    }

    let mut emission_buf = [0.0; Profile::MAX_ALPHABET_SIZE];
    let mut trans_buf = [0.0; Profile::NUM_STATE_TRANSITIONS - 1];
    let mut state = ParseState::Header;
    let mut model_state = ModelState::Comp;

    #[derive(Default)]
    pub struct ProfileBuilder {
        name: Option<String>,
        accession: Option<String>,
        composition: Option<[f32; Profile::MAX_ALPHABET_SIZE]>,
        mat_emissions: Vec<Option<[f32; Profile::MAX_ALPHABET_SIZE]>>,
        ins_emissions: Vec<Option<[f32; Profile::MAX_ALPHABET_SIZE]>>,
        transitions: Vec<Option<[f32; 7]>>,
        fwd_tau: Option<f32>,
        fwd_lambda: Option<f32>,
    }

    impl ProfileBuilder {
        pub fn name(&mut self, name: String) -> &mut Self {
            self.name = Some(name);
            self
        }

        pub fn accession(&mut self, accession: String) -> &mut Self {
            self.accession = Some(accession);
            self
        }

        pub fn length(&mut self, length: usize) -> &mut Self {
            self.mat_emissions.resize(length + 1, None);
            self.ins_emissions.resize(length + 1, None);
            self.transitions.resize(length + 1, None);
            self
        }

        pub fn mat_emission(&mut self, pos: usize, probs: [f32; 20]) -> &mut Self {
            self.mat_emissions[pos] = Some(probs);
            self
        }

        pub fn ins_emission(&mut self, pos: usize, probs: [f32; 20]) -> &mut Self {
            self.ins_emissions[pos] = Some(probs);
            self
        }

        pub fn transition(&mut self, pos: usize, probs: [f32; 7]) -> &mut Self {
            self.transitions[pos] = Some(probs);
            self
        }

        pub fn composition(&mut self, comp: [f32; 20]) -> &mut Self {
            self.composition = Some(comp);
            self
        }

        pub fn fwd_tau(&mut self, fwd_tau: f32) -> &mut Self {
            self.fwd_tau = Some(fwd_tau);
            self
        }

        pub fn fwd_lambda(&mut self, fwd_lambda: f32) -> &mut Self {
            self.fwd_lambda = Some(fwd_lambda);
            self
        }

        pub fn build(self) -> Profile {
            todo!()
        }
    }

    let mut pf = ProfileBuilder::default();
    let mut pos = 0;
    for line in bytes
        .split(|b| *b == b'\n')
        .filter(|b| !b.is_empty())
        .map(std::str::from_utf8)
    {
        let mut tokens = line?.splitn(2, char::is_whitespace);

        match state {
            ParseState::Header => {
                let flag = tokens.next().unwrap_or("").parse::<P7HeaderFlag>()?;
                let value = tokens
                    .next()
                    .ok_or(anyhow!("missing value for HMM flag: {}", flag.as_ref()))?;

                match flag {
                    P7HeaderFlag::Name => {
                        pf.name(value.to_string());
                    }
                    P7HeaderFlag::Accession => {
                        pf.accession(value.to_string());
                    }
                    P7HeaderFlag::Length => {
                        pf.length(
                            value.parse().with_context(|| {
                                format!("failed to parse HMM length: {}", value)
                            })?,
                        );
                    }
                    P7HeaderFlag::Alphabet => {
                        if value != "amino" {
                            bail!("non-amino HMM alphabet found: {}", value);
                        };
                    }
                    P7HeaderFlag::Stats => {
                        let stats: Vec<&str> = value.splitn(4, char::is_whitespace).collect();
                        if stats[1] == "FORWARD" {
                            pf.fwd_tau(stats[2].parse().with_context(|| {
                                format!("failed to parse HMM Forward tau: {}", value)
                            })?);
                            pf.fwd_lambda(stats[3].parse().with_context(|| {
                                format!("failed to parse HMM Forward lambda: {}", value)
                            })?);
                        }
                    }
                    P7HeaderFlag::Hmm => {
                        // let re = Regex::new(
                        //     r"^\s*A\s+C\s+D\s+E\s+F\s+G\s+H\s+I\s+K\s+L\s+M\s+N\s+P\s+Q\s+R\s+S\s+T\s+V\s+W\s+Y\s*$",
                        // );

                        state = ParseState::Model
                    }
                    _ => {
                        // ignore most flags for now
                    }
                }
            }
            ParseState::Model => {
                let (flag, values) = match tokens.next() {
                    Some(p) => (p, tokens.next().unwrap_or_default()),
                    None => continue,
                };

                match model_state {
                    ModelState::Comp => {
                        if flag == "COMPO" {
                            parse_floats_into::<{ Profile::MAX_ALPHABET_SIZE }>(
                                values,
                                &mut emission_buf,
                            )?;

                            pf.composition(emission_buf);
                            model_state = ModelState::Ins;
                        }
                    }
                    ModelState::Mat => {
                        let new_pos = flag.parse::<usize>()?;
                        if new_pos != pos + 1 {
                            bail!("non-linear HMM positions: {pos}->{new_pos}");
                        } else {
                            pos = new_pos
                        }

                        parse_floats_into::<{ Profile::MAX_ALPHABET_SIZE }>(
                            values,
                            &mut emission_buf,
                        )?;

                        pf.mat_emission(pos, emission_buf);
                        model_state = ModelState::Ins;
                    }
                    ModelState::Ins => {
                        parse_floats_into::<{ Profile::MAX_ALPHABET_SIZE }>(
                            values,
                            &mut emission_buf,
                        )?;

                        pf.ins_emission(pos, emission_buf);
                        model_state = ModelState::Trans;
                    }
                    ModelState::Trans => {
                        parse_floats_into::<{ Profile::NUM_STATE_TRANSITIONS - 1 }>(
                            values,
                            &mut trans_buf,
                        )?;

                        pf.transition(pos, trans_buf);
                        model_state = ModelState::Mat;
                    }
                }
            }
        }
    }

    Ok(pf.build())
}

fn sequence_from_fasta_record_bytes(bytes: &[u8]) -> anyhow::Result<Sequence> {
    let header_newline_pos = match bytes.iter().position(|&b| b == b'\n') {
        Some(pos) => pos,
        None => bail!("no newline in FASTA record"),
    };

    let header_bytes = &bytes[1..header_newline_pos];
    let sequence_bytes = &bytes[(header_newline_pos + 1)..];

    let name = String::from_utf8(
        header_bytes
            .iter()
            .take_while(|&&b| b != b' ')
            .cloned()
            .collect(),
    )?;

    let details = match name.len().cmp(&(header_newline_pos - 1)) {
        std::cmp::Ordering::Less => Some(String::from_utf8(
            header_bytes[(name.len() + 1)..]
                .iter()
                .take_while(|&&b| b != b'\n')
                .cloned()
                .collect(),
        )?),
        std::cmp::Ordering::Equal => None,
        std::cmp::Ordering::Greater => bail!("fasta name longer than header"),
    };

    let mut utf8_bytes = Vec::with_capacity(sequence_bytes.len());
    let mut digital_bytes = Vec::with_capacity(sequence_bytes.len());

    utf8_bytes.push(255);
    digital_bytes.push(255);

    sequence_bytes
        .iter()
        .filter(|&&b| b != b'\n')
        .try_for_each(|b| {
            utf8_bytes.push(*b);
            digital_bytes.push(match UTF8_TO_DIGITAL_AMINO.get(b) {
                Some(b) => *b,
                None => bail!("unknown byte"),
            });
            Ok(())
        })?;

    digital_bytes.push(Profile::NON_RESIDUE_IDX as u8);

    Ok(Sequence {
        name,
        details,
        length: utf8_bytes.len() - 1,
        digital_bytes,
        utf8_bytes,
    })
}

dyn_clone::clone_trait_object!(SequenceDatabase);
pub trait SequenceDatabase: dyn_clone::DynClone + Send + Sync {
    fn get(&mut self, name: &str) -> Option<Sequence>;
    fn len(&self) -> usize;
    fn iter(&self) -> SequenceDatabaseIter;
}

pub struct SequenceDatabaseIter<'a> {
    inner: Box<dyn SequenceDatabase>,
    names_iter: Box<dyn DoubleEndedIterator<Item = &'a str> + 'a>,
}

impl<'a> SequenceDatabaseIter<'a> {
    pub fn names(self) -> Vec<&'a str> {
        self.names_iter.collect()
    }
}

impl<'a> Iterator for SequenceDatabaseIter<'a> {
    type Item = Sequence;

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
impl<'a> DoubleEndedIterator for SequenceDatabaseIter<'a> {
    fn next_back(&mut self) -> Option<Self::Item> {
        self.names_iter.next_back().and_then(|n| self.inner.get(n))
    }
}

impl<'a> ExactSizeIterator for SequenceDatabaseIter<'a> {}

pub struct LexicalFastaIndex {
    offsets: IndexMap<String, FastaOffset>,
}

#[derive(Clone)]
pub struct FastaOffset {
    start: usize,
    name_len_bytes: usize,
    details_len_bytes: usize,
    seq_len_bytes: usize,
    seq_len: usize,
}

impl FastaOffset {
    pub fn new(start: usize) -> Self {
        Self {
            start,
            name_len_bytes: 0,
            details_len_bytes: 0,
            seq_len_bytes: 0,
            seq_len: 0,
        }
    }
}

impl LexicalFastaIndex {
    fn from_path<P: AsRef<Path>>(path: P, n: usize) -> anyhow::Result<Self> {
        let mut file = File::open(path.as_ref())?;
        let sz = file.metadata().unwrap().len();
        let chunk_sz = sz / n as u64;

        let mut block_starts: Vec<u64> = (0..n).map(|i| i as u64).map(|i| i * chunk_sz).collect();
        // 2^16 gives us a ~65KiB buffer
        let mut buffer = [0; 2 << 16];

        block_starts.iter_mut().for_each(|start| {
            file.seek(std::io::SeekFrom::Start(*start))
                .expect("failed to seek");

            let mut offset = 0;
            'outer: loop {
                match file.read(&mut buffer) {
                    Ok(bytes_read) => {
                        for b in buffer[0..bytes_read].iter() {
                            if *b == b'>' {
                                *start += offset;
                                break 'outer;
                            }
                            offset += 1;
                        }
                    }
                    _ => panic!("failed to read from buffer"),
                }
            }
        });

        let mut block_ends: Vec<u64> = block_starts.iter().skip(1).map(|b| b - 1).collect();
        block_ends.push(sz - 1);

        let mut files: Vec<_> = block_starts
            .into_iter()
            .zip(block_ends)
            .map(|(start, end)| {
                let mut f = File::open(path.as_ref())?;
                f.seek(std::io::SeekFrom::Start(start))?;
                Ok((f.take(end - start + 1), start))
            })
            .collect::<anyhow::Result<_>>()?;

        let mut indexes: Vec<_> = files
            .par_iter_mut()
            .map(|(file, start)| Self::new(file, Some(*start)))
            .collect();

        let mut index = indexes.remove(0);
        indexes.into_iter().for_each(|i| index.extend(i));

        Ok(index)
    }

    fn new<R: Read>(mut data: R, start: Option<u64>) -> Self {
        let start = start.unwrap_or(0) as usize;

        let mut offsets = IndexMap::new();
        enum ParseState {
            Name,
            Details,
            Seq,
            Process,
        }

        let mut buffer = [0];
        // simple check to make sure we are reading a fasta
        // TODO: need more checks for format validation
        //       during the entire parsing process
        if data.read(&mut buffer).unwrap() == 1 {
            assert!(buffer[0] == b'>')
        } else {
            panic!()
        }

        let mut state = ParseState::Name;
        // 2^20 gives us a ~1MiB buffer
        let mut buffer = vec![0; 2 << 20];
        let mut name = String::new();
        let mut total_bytes_read = 1 + start;
        let mut seq_line_cnt = 0;

        let mut offset = FastaOffset::new(start);

        while let Ok(bytes_read) = data.read(&mut buffer) {
            // when we read 0 bytes, its the EOF
            if bytes_read == 0 {
                offset.seq_len_bytes = total_bytes_read
                    - (offset.start + offset.name_len_bytes + offset.details_len_bytes);
                offset.seq_len = offset.seq_len_bytes - seq_line_cnt - 1;

                name.shrink_to_fit();
                offsets.insert(name.clone(), offset);
                break;
            }

            for (i, &byte) in buffer[..bytes_read].iter().enumerate() {
                let current_offset = total_bytes_read + i;
                match state {
                    ParseState::Name => match byte {
                        b' ' => {
                            offset.name_len_bytes = current_offset - offset.start;
                            state = ParseState::Details
                        }
                        b'\n' => {
                            offset.name_len_bytes = current_offset - offset.start;
                            offset.details_len_bytes = 0;
                            state = ParseState::Seq
                        }
                        _ => name.push(byte as char),
                    },
                    ParseState::Details => {
                        if byte == b'\n' {
                            offset.details_len_bytes =
                                current_offset - (offset.start + offset.name_len_bytes);
                            state = ParseState::Seq
                        }
                    }
                    ParseState::Seq => {
                        if byte == b'\n' {
                            seq_line_cnt += 1;
                        }

                        if byte == b'>' {
                            println!("{seq_line_cnt}");
                            offset.seq_len_bytes = current_offset
                                - (offset.start + offset.name_len_bytes + offset.details_len_bytes);
                            offset.seq_len = offset.seq_len_bytes - seq_line_cnt - 1;

                            state = ParseState::Process
                        }
                    }
                    ParseState::Process => {
                        name.shrink_to_fit();
                        offsets.insert(name.clone(), offset);
                        // -1 since we found the '>' on the previous byte
                        offset = FastaOffset::new(current_offset - 1);
                        name = String::new();
                        seq_line_cnt = 0;
                        name.push(byte as char);
                        state = ParseState::Name
                    }
                }
            }

            total_bytes_read += bytes_read;
        }

        offsets.shrink_to_fit();
        Self { offsets }
    }

    pub fn len(&self) -> usize {
        self.offsets.len()
    }

    fn offset_by_name(&self, name: &str) -> Option<FastaOffset> {
        self.offsets.get(name).cloned()
    }

    fn extend(&mut self, other: Self) {
        self.offsets.extend(other.offsets);
    }
}

pub struct Fasta {
    path: PathBuf,
    file: File,
    index: Arc<LexicalFastaIndex>,
    buffer: Vec<u8>,
}

impl Clone for Fasta {
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

impl Fasta {
    pub fn from_path_par<P: AsRef<Path>>(path: P, n: usize) -> anyhow::Result<Self> {
        let index = Arc::new(LexicalFastaIndex::from_path(path.as_ref(), n)?);
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
        let index = Arc::new(LexicalFastaIndex::new(&mut file, None));

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

    pub fn get(&mut self, name: &str) -> Option<Sequence> {
        let offset = self.index.offset_by_name(name)?;

        self.buffer.resize(
            offset.details_len_bytes + offset.name_len_bytes + offset.seq_len_bytes,
            0u8,
        );

        self.file
            .seek(std::io::SeekFrom::Start(offset.start as u64))
            .expect("failed to seek in Fasta::get()");

        self.file
            .read_exact(&mut self.buffer)
            .expect("failed to read in Fasta::get()");

        Some(
            sequence_from_fasta_record_bytes(&self.buffer).unwrap_or_else(|e| {
                panic!("failed to produce Sequence in Fasta::get()\nError: {e}");
            }),
        )
    }

    pub fn names_iter(&self) -> impl Iterator<Item = &str> {
        self.index.offsets.keys().map(|k| k.as_str())
    }

    pub fn lengths_iter(&self) -> impl Iterator<Item = usize> + use<'_> {
        self.index.offsets.values().map(|v| v.seq_len)
    }
}

impl SequenceDatabase for Fasta {
    fn get(&mut self, name: &str) -> Option<Sequence> {
        self.get(name)
    }

    fn len(&self) -> usize {
        self.len()
    }

    fn iter(&self) -> SequenceDatabaseIter {
        SequenceDatabaseIter {
            inner: Box::new(self.clone()),
            names_iter: Box::new(self.index.offsets.keys().map(|s| s.as_str())),
        }
    }
}
