mod rayon;

use std::{
    fs::File,
    io::{Read, Seek},
    path::{Path, PathBuf},
    sync::Arc,
};

use anyhow::bail;
use indexmap::IndexMap;
use libnail::{
    alphabet::UTF8_TO_DIGITAL_AMINO,
    structs::{Profile, Sequence},
};

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
    name_len: usize,
    details_len: usize,
    seq_len: usize,
}

impl FastaOffset {
    pub fn new(start: usize) -> Self {
        Self {
            start,
            name_len: 0,
            details_len: 0,
            seq_len: 0,
        }
    }
}

impl LexicalFastaIndex {
    fn new<R: Read>(mut data: R) -> Self {
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
        let mut buffer = [0; 8192];
        let mut name = String::new();
        let mut total_bytes_read = 1;

        let mut offset = FastaOffset::new(0);

        while let Ok(bytes_read) = data.read(&mut buffer) {
            // when we read 0 bytes, its the EOF
            if bytes_read == 0 {
                offset.seq_len =
                    total_bytes_read - (offset.start + offset.name_len + offset.details_len);
                name.shrink_to_fit();
                offsets.insert(name.clone(), offset);
                break;
            }

            for (i, &byte) in buffer[..bytes_read].iter().enumerate() {
                let current_offset = total_bytes_read + i;
                match state {
                    ParseState::Name => match byte {
                        b' ' => {
                            offset.name_len = current_offset - offset.start;
                            state = ParseState::Details
                        }
                        b'\n' => {
                            offset.name_len = current_offset - offset.start;
                            offset.details_len = 0;
                            state = ParseState::Seq
                        }
                        _ => name.push(byte as char),
                    },
                    ParseState::Details => {
                        if byte == b'\n' {
                            offset.details_len = current_offset - (offset.start + offset.name_len);
                            state = ParseState::Seq
                        }
                    }
                    ParseState::Seq => {
                        if byte == b'>' {
                            offset.seq_len = current_offset
                                - (offset.start + offset.name_len + offset.details_len);
                            state = ParseState::Process
                        }
                    }
                    ParseState::Process => {
                        name.shrink_to_fit();
                        offsets.insert(name.clone(), offset);
                        // -1 since we found the '>' on the previous byte
                        offset = FastaOffset::new(current_offset - 1);
                        name = String::new();
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
    pub fn from_path<P: AsRef<Path>>(path: P) -> anyhow::Result<Self> {
        let mut file = File::open(path.as_ref())?;
        let index = Arc::new(LexicalFastaIndex::new(&mut file));

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

        self.buffer
            .resize(offset.details_len + offset.name_len + offset.seq_len, 0u8);

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
