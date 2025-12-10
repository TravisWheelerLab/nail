use std::{
    fs::File,
    io::{Read, Seek, SeekFrom},
    path::{Path, PathBuf},
    sync::Arc,
};

use crate::io::{util::ByteBufferExt, DatabaseIter};

use libnail::structs::{
    profile::{ProfileBuilder, Transition},
    Profile,
};

use anyhow::{bail, Context};
use indexmap::IndexMap;
use strum::{AsRefStr, EnumIter, EnumString, IntoEnumIterator};

use super::{Database, DatabaseValues, Delimiter, Index, RecordParser};

#[derive(Clone, Default)]
pub struct P7HmmOffset {
    // points to the 'H' byte in the 'HMMER3/f ...' line
    start: usize,
    n_bytes: usize,
}

impl std::fmt::Debug for P7HmmOffset {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}->+{}", self.start, self.n_bytes)
    }
}

enum P7HmmParserState {
    Name,
    Delim,
}

pub struct P7HmmParser {
    state: P7HmmParserState,
    name: String,
    offset: P7HmmOffset,
}

impl RecordParser for P7HmmParser {
    const DELIM: &'static [u8] = b"//";
    const DELIM_TYPE: Delimiter = Delimiter::Terminating;
    type Offset = P7HmmOffset;
    type Record = Profile;

    fn new(start_pos: u64) -> Self {
        Self {
            state: P7HmmParserState::Name,
            name: String::new(),
            offset: P7HmmOffset {
                start: start_pos as usize,
                n_bytes: 0,
            },
        }
    }

    fn offset(&mut self, line: &[u8], line_start: u64) -> Option<(String, Self::Offset)> {
        use P7HmmParserState::*;

        self.offset.n_bytes += line.len();

        let first_word = match line.first_word() {
            Ok(w) => w,
            Err(_) => return None,
        };

        let mut ret = None;
        match (&self.state, first_word.parse()) {
            (Name, Ok(P7HeaderFlag::Name)) => {
                let name_start = line[P7HeaderFlag::Name.as_ref().len()..]
                    .iter()
                    .position(|b| !b.is_ascii_whitespace())
                    .map(|i| i + P7HeaderFlag::Name.as_ref().len())
                    .expect("missing name in p7HMM");

                self.name.push_str(
                    line.word_from(name_start)
                        .expect("failed to parse name in p7HMM line"),
                );

                self.state = Delim;
            }
            (Delim, Ok(P7HeaderFlag::Terminate)) => {
                ret = Some((self.name.as_str().to_owned(), self.offset.clone()));

                self.state = Name;
                self.name.clear();
                // if the line has a newline at the end, we want to add one
                let maybe_one = (line[line.len() - 1] != b'\n') as usize;
                self.offset.start = line_start as usize + line.len() + maybe_one;
                self.offset.n_bytes = 0;
            }
            _ => {}
        }
        ret
    }

    fn parse(buf: &[u8]) -> anyhow::Result<Self::Record> {
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
        let mut swap_buf = [0.0; Profile::NUM_STATE_TRANSITIONS - 1];
        let mut state = ParseState::Header;
        let mut model_state = ModelState::Comp;

        let mut prf = ProfileBuilder::default();
        let mut pos = 0;
        for line in buf.split(|b| *b == b'\n').filter(|b| !b.is_empty()) {
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
                                swap_buf[Transition::from(t) as usize] = trans_buf[t as usize]
                            });

                            prf.transition(pos, swap_buf);

                            model_state = ModelState::Mat;
                        }
                    }
                }
            }
        }

        prf.build()
    }
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
    #[strum(serialize = "//")]
    Terminate,
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

pub type P7HmmIndex = Index<IndexMap<String, P7HmmOffset>, P7HmmParser>;
pub struct P7Hmm {
    path: PathBuf,
    file: File,
    pub(crate) index: Arc<P7HmmIndex>,
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
    pub fn from_path<P: AsRef<Path>>(path: P) -> anyhow::Result<Self> {
        let index = Arc::new(P7HmmIndex::from_path::<20, 1>(path.as_ref())?);
        let file = File::open(path.as_ref())?;

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
        let offset = self.index.get(name)?;

        self.buffer.resize(offset.n_bytes, 0u8);

        self.file
            .seek(SeekFrom::Start(offset.start as u64))
            .expect("failed to seek in P7Hmm::get()");

        self.file
            .read_exact(&mut self.buffer)
            .expect("failed to read in P7Hmm::get()");

        Some(P7HmmParser::parse(&self.buffer).unwrap_or_else(|e| {
            panic!("failed to produce Profile in P7Hmm::get()\nError: {e}");
        }))
    }

    pub fn names_iter(&self) -> impl Iterator<Item = &str> {
        self.index.inner.keys().map(|k| k.as_str())
    }
}

impl Database<Profile> for P7Hmm {
    fn get(&mut self, name: &str) -> Option<Profile> {
        self.get(name)
    }

    fn len(&self) -> usize {
        self.len()
    }

    fn iter(&self) -> DatabaseIter<Profile> {
        DatabaseIter {
            inner: Box::new(self.clone()),
            names_iter: Box::new(self.index.inner.keys().map(|s| s.as_str())),
        }
    }

    fn values(&self) -> DatabaseValues<Profile> {
        DatabaseValues {
            inner: Box::new(self.clone()),
            names_iter: Box::new(self.index.inner.keys().map(|s| s.as_str())),
        }
    }
}

#[cfg(test)]
mod tests {
    use std::{fs::read_to_string, io::Cursor};

    use super::*;

    #[test]
    fn test_profile_serialize() -> anyhow::Result<()> {
        let mut bytes = vec![];
        let mut file = File::open(
            std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../fixtures/query.hmm"),
        )?;
        file.read_to_end(&mut bytes)?;
        let prf = P7HmmParser::parse(&bytes)?;

        let mut buf = vec![];

        prf.serialize(&mut buf)?;

        let de = Profile::deserialize(Cursor::new(buf))?;

        assert_eq!(prf.name, de.name);
        assert_eq!(prf.accession, de.accession);
        assert_eq!(prf.consensus_seq_bytes_utf8, de.consensus_seq_bytes_utf8);
        assert_eq!(prf.core_transitions, de.core_transitions);
        assert_eq!(prf.emission_scores[0], de.emission_scores[0]);
        Ok(())
    }

    #[test]
    fn test_p7hmm_index_starts() -> anyhow::Result<()> {
        let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../fixtures/query.hmm");
        let hmm_str = read_to_string(&path)?;
        let starts: Vec<usize> = hmm_str
            .as_bytes()
            .windows(10)
            .enumerate()
            .filter_map(|(i, w)| (w == b"HMMER3/f [").then_some(i))
            .collect();

        let names = hmm_str
            .lines()
            .filter(|l| l.starts_with(P7HeaderFlag::Name.as_ref()))
            .map(|l| l.split_whitespace().nth(1).unwrap())
            .collect::<Vec<&str>>();

        let index = P7HmmIndex::from_path::<15, 1>(&path)?;
        starts.iter().zip(names.iter()).for_each(|(s, n)| {
            let o = index.get(n).unwrap();
            assert_eq!(o.start, *s);
        });

        let index = P7HmmIndex::from_path::<15, 2>(&path)?;
        starts.iter().zip(names.iter()).for_each(|(s, n)| {
            let o = index.get(n).unwrap();
            assert_eq!(o.start, *s);
        });

        let index = P7HmmIndex::from_path::<15, 3>(&path)?;
        starts.iter().zip(names.iter()).for_each(|(s, n)| {
            let o = index.get(n).unwrap();
            assert_eq!(o.start, *s);
        });

        let index = P7HmmIndex::from_path::<15, 4>(&path)?;
        starts.iter().zip(names.iter()).for_each(|(s, n)| {
            let o = index.get(n).unwrap();
            assert_eq!(o.start, *s);
        });

        let index = P7HmmIndex::from_path::<15, 20>(&path)?;
        starts.iter().zip(names.iter()).for_each(|(s, n)| {
            let o = index.get(n).unwrap();
            assert_eq!(o.start, *s);
        });

        Ok(())
    }
}
