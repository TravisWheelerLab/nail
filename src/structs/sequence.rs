use seq_io::fasta::{Reader, Record};
use std::path::Path;

use crate::alphabet::UTF8_TO_DIGITAL_AMINO;
use anyhow::{Context, Result};
use thiserror::Error;

#[derive(Error, Debug)]
#[error("unknown fasta sequence character")]
struct UnknownSequenceCharacterError;

/// This holds the both the "digital" data and string data of a biological sequence.
pub struct Sequence {
    /// The name of the sequence
    pub name: String,
    /// The length of the sequence
    pub length: usize,
    /// The "digital" data of the sequence. These are the string bytes, but mapped to [0u8..25u8]
    pub digital_bytes: Vec<u8>,
    /// The string data of the sequence. These are the UTF8 bytes that make up the sequence in the "normal" alphabet
    pub utf8_bytes: Vec<u8>,
}

impl Sequence {
    pub fn amino_from_fasta<P: AsRef<Path>>(path: P) -> Result<Vec<Self>> {
        let mut seqs: Vec<Self> = vec![];

        let mut reader = Reader::from_path(path).unwrap();

        while let Some(record) = reader.next() {
            let record = record.expect("Error reading record");
            let record_name = String::from_utf8(record.head().to_vec())?;
            // We want position 1 of the sequence to be at index 1, so we'll buffer with 255
            let mut utf8_bytes: Vec<u8> = vec![255];
            let mut digital_bytes: Vec<u8> = vec![255];

            for line in record.seq_lines() {
                for utf8_byte in line {
                    utf8_bytes.push(*utf8_byte);

                    let amino_byte = match UTF8_TO_DIGITAL_AMINO.get(utf8_byte) {
                        Some(b) => b,
                        None => {
                            return Err(UnknownSequenceCharacterError).with_context(|| {
                                format!("failed to parse utf8 byte: {}", utf8_byte)
                            })
                        }
                    };
                    digital_bytes.push(*amino_byte)
                }
            }

            seqs.push(Sequence {
                name: record_name,
                length: digital_bytes.len() - 1,
                digital_bytes,
                utf8_bytes,
            });
        }
        Ok(seqs)
    }
}
