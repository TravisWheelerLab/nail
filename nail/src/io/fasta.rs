use std::{
    fs::File,
    io::{BufWriter, Read, Seek, Write},
    path::{Path, PathBuf},
    sync::{Arc, Mutex},
};

use anyhow::{anyhow, bail};
use indexmap::IndexMap;
use libnail::{
    alphabet::UTF8_TO_DIGITAL_AMINO,
    structs::{Profile, Sequence},
    util::IterPrint,
};
use rayon::iter::{IntoParallelRefMutIterator, ParallelIterator};

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
    pub(crate) inner: Box<dyn SequenceDatabase>,
    pub(crate) names_iter: Box<dyn DoubleEndedIterator<Item = &'a str> + 'a>,
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

#[derive(Default)]
pub struct LexicalFastaIndex {
    pub(crate) offsets: IndexMap<String, FastaOffset>,
}

#[derive(Clone)]
pub struct FastaOffset {
    // this points to the '>' byte that starts the fasta record
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
            'inner: loop {
                match file.read(&mut buffer) {
                    Ok(0) => {
                        *start += offset;
                        break 'inner;
                    }
                    Ok(bytes_read) => {
                        for b in buffer[0..bytes_read].iter() {
                            if *b == b'>' {
                                *start += offset;
                                break 'inner;
                            }
                            offset += 1;
                        }
                    }
                    _ => panic!("failed to read from buffer"),
                }
            }
        });

        // if we have duplicate starts, then the file was
        // really small and we had at least one thread
        // reach EOF before ever finding it's own chunk
        block_starts.dedup();

        match block_starts.last() {
            // if the last start points to EOF,
            // then we don't want it around
            Some(st) if *st == sz => {
                block_starts.pop();
            }
            Some(_) => {}
            None => bail!("empty block_starts"),
        }

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

        let index = Arc::new(Mutex::new(Self::default()));

        files
            .par_iter_mut()
            .try_for_each(|(file, start)| Self::build(file, index.clone(), Some(*start)))?;

        Ok(Arc::try_unwrap(index)
            .ok()
            .expect("other Arc clones exist")
            .into_inner()
            .unwrap())
    }

    fn build<R: Read>(
        mut data: R,
        index: Arc<Mutex<Self>>,
        start: Option<u64>,
    ) -> anyhow::Result<()> {
        const N: usize = 1_000;
        let start = start.unwrap_or(0) as usize;

        enum ParseState {
            Name,
            Details,
            Seq,
            Process,
        }

        let mut buf = [0];
        // simple check to make sure we are reading a fasta
        // TODO: need more checks for format validation
        //       during the entire parsing process
        if data.read(&mut buf).unwrap() == 1 {
            assert!(buf[0] == b'>')
        } else {
            panic!()
        }

        let mut offsets: Vec<(String, FastaOffset)> = Vec::with_capacity(N);

        let mut state = ParseState::Name;
        // 2^20 gives us a ~1MiB buffer
        let mut buf = vec![0; 2 << 20];
        let mut name = String::new();
        let mut total_bytes_read = 1 + start;
        let mut seq_line_cnt = 0;

        let mut offset = FastaOffset::new(start);

        let process_fn = |mut name: String,
                          offset: FastaOffset,
                          offsets: &mut Vec<(String, FastaOffset)>,
                          done: bool| {
            if name.trim().is_empty() {
                bail!("failed to parse name from fasta buffer");
            }

            name.shrink_to_fit();
            offsets.push((name, offset));

            match index.try_lock() {
                Ok(mut guard) => {
                    guard.offsets.extend(offsets.drain(..));
                }
                Err(std::sync::TryLockError::WouldBlock) => {
                    if offsets.len() == N || done {
                        index
                            .lock()
                            .map_err(|_| anyhow!("mutex poisoned"))?
                            .offsets
                            .extend(offsets.drain(..));
                    }
                }
                Err(_) => bail!("mutex poisoned"),
            }

            Ok(())
        };

        while let Ok(bytes_read) = data.read(&mut buf) {
            // TODO (IMPORTANT): truncate to final newline & seek back

            // when we read 0 bytes, its the EOF
            if bytes_read == 0 {
                offset.seq_len_bytes = total_bytes_read
                    - (offset.start + offset.name_len_bytes + offset.details_len_bytes);
                offset.seq_len = offset.seq_len_bytes - seq_line_cnt - 1;

                process_fn(name, offset, &mut offsets, true)?;
                break;
            }

            for (i, &byte) in buf[..bytes_read].iter().enumerate() {
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
                            offset.seq_len_bytes = current_offset
                                - (offset.start + offset.name_len_bytes + offset.details_len_bytes);
                            offset.seq_len = offset.seq_len_bytes - seq_line_cnt - 1;

                            state = ParseState::Process
                        }
                    }
                    ParseState::Process => {
                        process_fn(name, offset, &mut offsets, false)?;
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
        Ok(())
    }

    pub fn len(&self) -> usize {
        self.offsets.len()
    }

    fn offset_by_name(&self, name: &str) -> Option<FastaOffset> {
        self.offsets.get(name).cloned()
    }
}

pub struct Fasta {
    path: PathBuf,
    file: File,
    pub(crate) index: Arc<LexicalFastaIndex>,
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

    pub fn len(&self) -> usize {
        self.index.len()
    }

    pub fn write<W: Write>(&self, out: W) -> anyhow::Result<()> {
        let mut out = BufWriter::new(out);
        self.iter().try_for_each(|s| writeln!(out, "{s}"))?;
        Ok(())
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
