use seq_io::fasta::{Reader, Record};
use std::fmt::{Debug, Formatter};
use std::path::Path;

use crate::alphabet::{AMINO_INVERSE_MAP, UTF8_TO_DIGITAL_AMINO};
use anyhow::{Context, Result};
use thiserror::Error;

#[derive(Error, Debug)]
#[error("unknown UTF8 sequence byte: {byte}")]
pub struct UnknownUtf8SequenceByteError {
    byte: u8,
}

#[derive(Error, Debug)]
#[error("unknown digital sequence byte: {byte}")]
pub struct UnknownDigitalSequenceByteError {
    byte: u8,
}

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
            let record_header = String::from_utf8(record.head().to_vec())?;
            let record_name = record_header.split_whitespace().next().unwrap().to_string();
            // We want position 1 of the sequence to be at index 1, so we'll buffer with 255
            let mut utf8_bytes: Vec<u8> = vec![255];
            let mut digital_bytes: Vec<u8> = vec![255];

            for line in record.seq_lines() {
                for utf8_byte in line {
                    utf8_bytes.push(*utf8_byte);

                    let digital_byte = match UTF8_TO_DIGITAL_AMINO.get(utf8_byte) {
                        Some(b) => b,
                        None => {
                            return Err(UnknownUtf8SequenceByteError { byte: *utf8_byte }.into())
                        }
                    };
                    digital_bytes.push(*digital_byte)
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

    pub fn from_digital(bytes: &[u8]) -> Result<Self> {
        let mut digital_bytes: Vec<u8> = vec![255; bytes.len() + 1];
        digital_bytes[1..].copy_from_slice(bytes);
        let mut utf8_bytes: Vec<u8> = vec![255; digital_bytes.len()];

        for (idx, digital_byte) in digital_bytes[1..].iter().enumerate() {
            let utf8_byte = match AMINO_INVERSE_MAP.get(digital_byte) {
                Some(b) => *b,
                None => {
                    return Err(UnknownDigitalSequenceByteError {
                        byte: *digital_byte,
                    }
                    .into())
                }
            };
            utf8_bytes[idx + 1] = utf8_byte;
        }

        Ok(Sequence {
            name: "".to_string(),
            length: utf8_bytes.len() - 1,
            digital_bytes,
            utf8_bytes,
        })
    }

    pub fn from_utf8(bytes: &[u8]) -> Result<Self> {
        let mut utf8_bytes: Vec<u8> = vec![255; bytes.len() + 1];
        utf8_bytes[1..].copy_from_slice(bytes);
        let mut digital_bytes: Vec<u8> = vec![255; utf8_bytes.len()];

        for (idx, utf8_byte) in utf8_bytes[1..].iter().enumerate() {
            let digital_byte = match UTF8_TO_DIGITAL_AMINO.get(utf8_byte) {
                Some(b) => *b,
                None => return Err(UnknownUtf8SequenceByteError { byte: *utf8_byte }.into()),
            };
            digital_bytes[idx + 1] = digital_byte;
        }

        Ok(Sequence {
            name: "".to_string(),
            length: digital_bytes.len() - 1,
            digital_bytes,
            utf8_bytes,
        })
    }
}

impl Debug for Sequence {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", std::str::from_utf8(&self.utf8_bytes).unwrap())?;
        Ok(())
    }
}
