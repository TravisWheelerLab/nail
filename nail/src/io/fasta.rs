use anyhow::bail;
use libnail::{
    alphabet::UTF8_TO_DIGITAL_AMINO,
    structs::{Profile, Sequence},
};

use crate::io::{ByteBufferExt, Offset};

use super::{Delimiter, RecordParser};

#[derive(Clone)]
pub struct FastaParser {
    name: String,
    offset: Offset,
}

impl RecordParser for FastaParser {
    const DELIM: &'static [u8] = b">";
    const DELIM_TYPE: Delimiter = Delimiter::Initiating;
    type Record = Sequence;

    fn new(_: u64) -> Self {
        Self {
            name: String::new(),
            offset: Offset::default(),
        }
    }

    fn offset(&mut self, line: &[u8], line_start: u64) -> Option<(String, Offset)> {
        let mut ret = None;

        if line[0] == Self::DELIM[0] {
            // this means we've hit a new record, so we're going to:
            //   1. return the last offset we built, and
            //   2. start building the next offset

            ret = Some((self.name.as_str().to_owned(), self.offset.clone()));

            self.name.clear();
            self.name
                .push_str((&line[1..]).first_word().unwrap_or_default());
            self.offset.start = line_start;
            self.offset.n_bytes = 0;
        }

        self.offset.n_bytes += line.len();

        ret
    }

    fn parse(buf: &[u8]) -> anyhow::Result<Self::Record> {
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
}

// #[cfg(test)]
// mod tests {
//     use std::fs::read_to_string;

//     use super::*;

//     #[test]
//     fn test_fasta_index_starts() -> anyhow::Result<()> {
//         let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../fixtures/target.fa");
//         let fasta_str = read_to_string(&path)?;
//         let starts: Vec<usize> = fasta_str
//             .as_bytes()
//             .iter()
//             .enumerate()
//             .filter_map(|(i, b)| (*b == b'>').then_some(i))
//             .collect();

//         let fasta_bytes = fasta_str.as_bytes();
//         let names: Vec<&str> = starts
//             .iter()
//             .filter_map(|p| fasta_bytes.word_from(p + 1).ok())
//             .collect();

//         let index = FastaIndex::from_path::<15, 1>(&path)?;
//         starts.iter().zip(names.iter()).for_each(|(s, n)| {
//             let o = index.get(n).unwrap();
//             assert_eq!(o.start, *s as u64);
//         });

//         let index = FastaIndex::from_path::<15, 2>(&path)?;
//         starts.iter().zip(names.iter()).for_each(|(s, n)| {
//             let o = index.get(n).unwrap();
//             assert_eq!(o.start, *s as u64);
//         });

//         let index = FastaIndex::from_path::<15, 3>(&path)?;
//         starts.iter().zip(names.iter()).for_each(|(s, n)| {
//             let o = index.get(n).unwrap();
//             assert_eq!(o.start, *s as u64);
//         });

//         let index = FastaIndex::from_path::<15, 4>(&path)?;
//         starts.iter().zip(names.iter()).for_each(|(s, n)| {
//             let o = index.get(n).unwrap();
//             assert_eq!(o.start, *s as u64);
//         });

//         let index = FastaIndex::from_path::<15, 20>(&path)?;
//         starts.iter().zip(names.iter()).for_each(|(s, n)| {
//             let o = index.get(n).unwrap();
//             assert_eq!(o.start, *s as u64);
//         });

//         Ok(())
//     }
// }
