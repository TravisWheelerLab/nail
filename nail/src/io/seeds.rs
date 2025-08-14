use std::{
    fs::File,
    io::{Read, Seek, SeekFrom, Write},
    path::{Path, PathBuf},
    sync::Arc,
};

use super::{Database, DatabaseIter, Delimiter, Index, RecordParser};
use crate::io::{ByteBufferExt, IndexInner};

use anyhow::bail;
use indexmap::IndexMap;
use libnail::align::structs::Seed;

#[derive(Default, Clone, Debug, PartialEq)]
pub struct SeedOffset {
    start: u64,
    n_bytes: usize,
}

pub struct SeedParser;

impl RecordParser for SeedParser {
    const DELIM: &'static [u8] = b"\n";
    const DELIM_TYPE: Delimiter = Delimiter::Terminating;
    type Offset = SeedOffset;
    type Record = SeedList;

    fn new(_: u64) -> Self {
        Self
    }

    fn offset(&mut self, line: &[u8], line_start: u64) -> Option<(String, Self::Offset)> {
        let name = line.first_word().ok()?.to_string();

        Some((
            name,
            SeedOffset {
                start: line_start,
                n_bytes: line.len(),
            },
        ))
    }

    fn parse(buf: &[u8]) -> anyhow::Result<Self::Record> {
        let mut seeds = vec![];

        let buf = match buf.last() {
            Some(b'\n') => &buf[..buf.len() - 1],
            Some(_) => buf,
            None => bail!("empty buffer"),
        };

        for line in buf.split(|b| *b == b'\n') {
            let line = line.as_str()?;
            let tokens = line.split_whitespace().collect::<Vec<_>>();
            seeds.push((
                tokens[1].into(),
                Seed {
                    seq_start: tokens[4].parse()?,
                    seq_end: tokens[5].parse()?,
                    prf_start: tokens[2].parse()?,
                    prf_end: tokens[3].parse()?,
                    score: tokens[6].parse()?,
                },
            ))
        }
        Ok(seeds)
    }
}

pub struct SeedsIndexInner {
    index: IndexMap<String, SeedOffset>,
}

impl IndexInner<SeedOffset> for SeedsIndexInner {
    fn new() -> Self {
        Self {
            index: IndexMap::new(),
        }
    }

    fn len(&self) -> usize {
        self.index.len()
    }

    fn get(&self, name: &str) -> Option<&SeedOffset> {
        self.index.get(name)
    }

    fn extend<T>(&mut self, iter: T)
    where
        T: IntoIterator<Item = (String, SeedOffset)>,
    {
        iter.into_iter().for_each(|(name, offset)| {
            self.index
                .entry(name)
                .and_modify(|existing| {
                    let (s1, e1) = (existing.start, existing.start + existing.n_bytes as u64);
                    let (s2, e2) = (offset.start, offset.start + offset.n_bytes as u64);

                    existing.start = s1.min(s2);
                    existing.n_bytes = (e1.max(e2) - existing.start) as usize;
                })
                .or_insert(offset);
        });
    }
}

pub type SeedsIndex = Index<SeedsIndexInner, SeedParser>;
pub type SeedList = Vec<(String, Seed)>;

pub struct Seeds {
    path: PathBuf,
    file: File,
    pub(crate) index: Arc<SeedsIndex>,
    buffer: Vec<u8>,
}

impl Clone for Seeds {
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

impl Seeds {
    pub fn from_path<P: AsRef<Path>>(path: P) -> anyhow::Result<Self> {
        let index = Arc::new(Index::from_path::<15, 4>(path.as_ref())?);
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

    pub fn write<W: Write>(&self, out: W) -> anyhow::Result<()> {
        todo!()
    }

    pub fn get(&mut self, name: &str) -> Option<SeedList> {
        let offset = self.index.get(name)?;

        self.buffer.resize(offset.n_bytes, 0u8);

        self.file
            .seek(SeekFrom::Start(offset.start))
            .expect("failed to seek in Seeds::get()");

        self.file
            .read_exact(&mut self.buffer)
            .expect("failed to read in Seeds::get()");

        Some(SeedParser::parse(&self.buffer).unwrap_or_else(|e| {
            panic!("failed to produce Seeds in Seeds::get()\nError: {e}");
        }))
    }

    pub fn names_iter(&self) -> impl Iterator<Item = &str> {
        self.index.inner.index.keys().map(|k| k.as_str())
    }
}

impl Database<SeedList> for Seeds {
    fn get(&mut self, name: &str) -> Option<Vec<(String, Seed)>> {
        self.get(name)
    }

    fn len(&self) -> usize {
        self.len()
    }

    fn iter(&self) -> DatabaseIter<SeedList> {
        DatabaseIter {
            inner: Box::new(self.clone()),
            names_iter: Box::new(self.index.inner.index.keys().map(|s| s.as_str())),
        }
    }
}

#[cfg(test)]
mod tests {
    use std::{collections::HashMap, fs::read_to_string};

    use super::*;

    #[test]
    fn test_seeds_index() -> anyhow::Result<()> {
        let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../fixtures/a.seeds");
        let seeds_str = read_to_string(&path)?;

        let offsets = {
            let mut starts = HashMap::new();
            let mut ends = HashMap::new();

            let mut pos = 0usize;
            seeds_str.lines().for_each(|line| {
                let name = line.split_whitespace().next().unwrap();
                starts.entry(name).or_insert(pos);
                pos += line.len() + 1;
                ends.entry(name)
                    .and_modify(|n: &mut usize| *n = (*n).max(pos))
                    .or_insert(pos);
            });

            starts
                .iter()
                .map(|(n, s)| {
                    let e = ends.get(n).unwrap();
                    (
                        n.to_string(),
                        SeedOffset {
                            start: *s as u64,
                            n_bytes: e - s,
                        },
                    )
                })
                .collect::<Vec<_>>()
        };

        let index = SeedsIndex::from_path::<15, 1>(&path)?;
        offsets.iter().for_each(|(n, o)| {
            let o2 = index.get(n).unwrap();
            assert_eq!(o, o2);
        });

        let index = SeedsIndex::from_path::<15, 2>(&path)?;
        offsets.iter().for_each(|(n, o)| {
            let o2 = index.get(n).unwrap();
            assert_eq!(o, o2);
        });

        let index = SeedsIndex::from_path::<15, 3>(&path)?;
        offsets.iter().for_each(|(n, o)| {
            let o2 = index.get(n).unwrap();
            assert_eq!(o, o2);
        });

        let index = SeedsIndex::from_path::<15, 4>(&path)?;
        offsets.iter().for_each(|(n, o)| {
            let o2 = index.get(n).unwrap();
            assert_eq!(o, o2);
        });

        let index = SeedsIndex::from_path::<15, 20>(&path)?;
        offsets.iter().for_each(|(n, o)| {
            let o2 = index.get(n).unwrap();
            assert_eq!(o, o2);
        });

        Ok(())
    }
}
