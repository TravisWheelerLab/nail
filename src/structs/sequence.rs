use crate::alphabet::AMINO_MAP;

use seq_io::fasta::Reader;
use std::path::Path;

use anyhow::{Context, Result};
use thiserror::Error;

#[derive(Error, Debug)]
#[error("unknown fasta sequence character")]
struct UnknownSequenceCharacterError;

pub struct Sequence {
    pub length: usize,
    pub data: Vec<u8>,
}

impl Sequence {
    pub fn from_fasta<P: AsRef<Path>>(path: P) -> Result<Vec<Self>> {
        let mut seqs: Vec<Self> = vec![];

        let mut reader = Reader::from_path(path).unwrap();

        while let Some(record) = reader.next() {
            let record = record.expect("Error reading record");
            // I think we want position 1 of the sequence to be at index 1, so we'll buffer with 255
            let mut amino_bytes: Vec<u8> = vec![255];
            for line in record.seq_lines() {
                for utf8_byte in line {
                    let amino_byte = match AMINO_MAP.get(&(*utf8_byte as char)) {
                        Some(b) => b,
                        None => {
                            return Err(UnknownSequenceCharacterError)
                                .with_context(|| format!("failed to parse byte: {}", utf8_byte))
                        }
                    };
                    amino_bytes.push(*amino_byte)
                }
            }

            seqs.push(Sequence {
                length: amino_bytes.len() - 1,
                data: amino_bytes,
            });
        }
        Ok(seqs)
    }
}
