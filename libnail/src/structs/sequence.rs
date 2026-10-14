use rand::Rng;
use rand_pcg::Lcg128Xsl64;
use std::fmt::{Debug, Display, Formatter};

use crate::alphabet::{
    AMINO_BACKGROUND_FREQUENCIES, AMINO_INVERSE_MAP, AMINO_SENTINEL, UTF8_TO_DIGITAL_AMINO,
};
use anyhow::{bail, Result};
use thiserror::Error;

use super::Profile;

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
#[derive(Default, PartialEq)]
pub struct Sequence {
    /// The name of the sequence
    pub name: String,
    /// The sequence details. If the sequence comes from a fasta, this
    /// is the information following the sequence name in the header
    pub details: Option<String>,
    /// The length of the sequence
    pub length: usize,
    /// The "digital" data of the sequence. These are
    /// the string bytes, but mapped to [0u8..25u8]
    pub digital_bytes: Vec<u8>,
    /// The string data of the sequence. These are the UTF8 bytes
    /// that make up the sequence in the "normal" alphabet
    pub utf8_bytes: Vec<u8>,
}

impl Sequence {
    pub fn random_amino(length: usize, rng: &mut Lcg128Xsl64) -> Self {
        let mut digital_bytes = (0..=length)
            .map(|_| {
                let roll: f32 = rng.gen();

                let mut sum = 0.0f32;
                let mut choice = 0u8;
                for (residue_idx, p) in AMINO_BACKGROUND_FREQUENCIES.iter().enumerate() {
                    sum += p;
                    if roll <= sum {
                        choice = residue_idx as u8;
                        break;
                    }
                }
                choice
            })
            .collect::<Vec<u8>>();

        // we should have <length + 1> bytes, since
        // we start with the sentinel character
        debug_assert_eq!(length + 1, digital_bytes.len());

        digital_bytes[0] = AMINO_SENTINEL;

        let mut utf8_bytes: Vec<u8> = digital_bytes
            .iter()
            .map(|b| *AMINO_INVERSE_MAP.get(b).expect(""))
            .collect();

        digital_bytes.push(Profile::NON_RESIDUE_IDX as u8);

        digital_bytes.shrink_to_fit();
        utf8_bytes.shrink_to_fit();

        Self {
            name: "random-seq".to_string(),
            details: None,
            length,
            digital_bytes,
            utf8_bytes,
        }
    }

    /// Build a sequence from the bytes of one fasta record.
    pub fn from_fasta(buf: &[u8]) -> anyhow::Result<Sequence> {
        let header_newline_pos = match buf.iter().position(|&b| b == b'\n') {
            Some(pos) => pos,
            None => bail!("no newline in FASTA record"),
        };

        let header_bytes = &buf[1..header_newline_pos];
        let sequence_bytes = &buf[(header_newline_pos + 1)..];

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
                    None => {
                        bail!("byte: \"{b}\" is not in the compressed amino acid alphabet")
                    }
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
            details: None,
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
            details: None,
            length: digital_bytes.len() - 1,
            digital_bytes,
            utf8_bytes,
        })
    }
}

impl Display for Sequence {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(f, ">{}", self.name)?;

        if let Some(ref details) = self.details {
            write!(f, " {details}")?
        };

        writeln!(f)?;

        // note: the utf8 bytes start with a padding byte of 255
        let mut iter = self.utf8_bytes[1..].chunks(80).peekable();

        while let Some(byte_chunk) = iter.next() {
            match std::str::from_utf8(byte_chunk) {
                Ok(seq_line) => {
                    write!(f, "{}", seq_line)?;
                    if iter.peek().is_some() {
                        // if we're not on the last
                        // line, add a linebreak
                        writeln!(f)?;
                    }
                }
                Err(_) => return Err(std::fmt::Error),
            }
        }
        Ok(())
    }
}

impl Debug for Sequence {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", std::str::from_utf8(&self.utf8_bytes[1..]).unwrap())?;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use rand::SeedableRng;
    use rand_pcg::Pcg64;

    use crate::{alphabet::AMINO_BACKGROUND_FREQUENCIES, structs::Sequence};

    #[test]
    fn test_random_amino() -> anyhow::Result<()> {
        const TOLERANCE: f32 = 1e-3;
        let mut counts = [0usize; 20];

        let mut rng = Pcg64::seed_from_u64(0);

        (0..1_000).for_each(|_| {
            let s = Sequence::random_amino(100, &mut rng);
            s.digital_bytes
                .iter()
                .skip(1)
                .take(100)
                .for_each(|&b| counts[b as usize] += 1);
        });

        let total_residues = counts.iter().sum::<usize>() as f32;

        counts
            .iter()
            .map(|&c| c as f32 / total_residues)
            .zip(AMINO_BACKGROUND_FREQUENCIES)
            .for_each(|(computed_f, correct_f)| {
                // println!("{correct_f:.3} {computed_f:.3}");
                let diff = (computed_f - correct_f).abs();
                assert!(
                    diff <= TOLERANCE,
                    "difference in computed residue frequencies is above tolerance: {diff} > {TOLERANCE}",
                );
            });
        Ok(())
    }
}
