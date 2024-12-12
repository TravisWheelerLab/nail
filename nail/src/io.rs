use std::{
    collections::HashMap,
    fs::File,
    io::{Read, Seek},
    path::{Path, PathBuf},
    sync::Arc,
};

use anyhow::bail;
use libnail::{alphabet::UTF8_TO_DIGITAL_AMINO, structs::Sequence};

pub trait SequenceDatabase: dyn_clone::DynClone + Send + Sync {
    fn get(&mut self, name: &str) -> anyhow::Result<Sequence>;
}

dyn_clone::clone_trait_object!(SequenceDatabase);

fn sequence_from_fasta_record_bytes(bytes: &[u8]) -> anyhow::Result<Sequence> {
    let header_end = match bytes.iter().position(|&b| b == b'\n') {
        Some(pos) => pos,
        None => bail!("no newline in FASTA record"),
    };

    let header_bytes = &bytes[1..header_end];
    let sequence_bytes = &bytes[(header_end + 1)..];

    let name = String::from_utf8(
        header_bytes
            .iter()
            .take_while(|&&b| b != b' ')
            .cloned()
            .collect(),
    )?;

    let details = match name.len().cmp(&header_end) {
        std::cmp::Ordering::Less => Some(String::from_utf8(
            header_bytes[(name.len() + 1)..]
                .iter()
                .take_while(|&&b| b != b'\n')
                .cloned()
                .collect(),
        )?),
        std::cmp::Ordering::Equal => None,
        std::cmp::Ordering::Greater => bail!("unexpected"),
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

    Ok(Sequence {
        name,
        details,
        length: utf8_bytes.len() - 1,
        digital_bytes,
        utf8_bytes,
    })
}

pub struct LexicalFastaIndex {
    offsets: HashMap<String, FastaOffset>,
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
        let mut offsets = HashMap::new();
        let mut buffer = [0; 8192];
        let mut name = String::new();
        let mut total_bytes_read = 1;

        let mut offset = FastaOffset::new(0);

        while let Ok(bytes_read) = data.read(&mut buffer) {
            // when we read 0 bytes, its the EOF
            if bytes_read == 0 {
                offset.seq_len =
                    total_bytes_read - (offset.start + offset.name_len + offset.details_len);
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
                        offsets.insert(name.clone(), offset);
                        // -1 since we found the '>' on the previous byte
                        offset = FastaOffset::new(current_offset - 1);
                        name.clear();
                        name.push(byte as char);
                        state = ParseState::Name
                    }
                }
            }

            total_bytes_read += bytes_read;
        }

        Self { offsets }
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

    pub fn get(&mut self, name: &str) -> anyhow::Result<Sequence> {
        let offset = self.index.offset_by_name(name).unwrap();

        self.buffer
            .resize(offset.details_len + offset.name_len + offset.seq_len, 0u8);
        self.file
            .seek(std::io::SeekFrom::Start(offset.start as u64))?;
        self.file.read_exact(&mut self.buffer)?;

        sequence_from_fasta_record_bytes(&self.buffer)
    }
}

impl SequenceDatabase for Fasta {
    fn get(&mut self, name: &str) -> anyhow::Result<Sequence> {
        self.get(name)
    }
}
