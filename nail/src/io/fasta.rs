use std::path::Path;

use libnail::structs::Sequence;
use libsail::format::Format;
use rayon::prelude::*;

use crate::io::Records;

/// A fasta file and where each of its records sits in it.
pub struct Fasta(Records);

impl Fasta {
    pub fn from_path(path: impl AsRef<Path>) -> anyhow::Result<Self> {
        Ok(Self(Records::map(path, Format::Fasta)?))
    }

    pub fn len(&self) -> usize {
        self.0.len()
    }

    /// The name of the `n`th record.
    pub fn name(&self, n: usize) -> anyhow::Result<String> {
        self.0.name(n)
    }

    /// The `n`th record.
    pub fn get(&self, n: usize) -> anyhow::Result<Sequence> {
        self.0.parse(n, Sequence::from_fasta)
    }

    pub fn iter(&self) -> impl ExactSizeIterator<Item = anyhow::Result<Sequence>> + '_ {
        (0..self.len()).map(|n| self.get(n))
    }

    pub fn par_iter(&self) -> impl IndexedParallelIterator<Item = anyhow::Result<Sequence>> + '_ {
        (0..self.len()).into_par_iter().map(|n| self.get(n))
    }
}

#[cfg(test)]
pub mod tests {
    use std::{fs::read_to_string, path::Path};

    use rayon::iter::ParallelIterator;

    use super::*;

    #[test]
    fn test_fasta_iter_parse() -> anyhow::Result<()> {
        let path = Path::new(env!("CARGO_MANIFEST_DIR")).join("../fixtures/target.fa");
        let database = Fasta::from_path(&path)?;
        let _ = database.iter().collect::<anyhow::Result<Vec<_>>>()?;
        Ok(())
    }

    #[test]
    fn test_fasta_records_in_file_order() -> anyhow::Result<()> {
        let path = Path::new(env!("CARGO_MANIFEST_DIR")).join("../fixtures/target.fa");

        // the name of each record, by a plain scan of the header lines
        let names: Vec<String> = read_to_string(&path)?
            .lines()
            .filter_map(|line| line.strip_prefix('>'))
            .map(|header| header.split(' ').next().unwrap().to_string())
            .collect();

        let fasta = Fasta::from_path(&path)?;
        assert_eq!(fasta.len(), names.len());

        for (n, name) in names.iter().enumerate() {
            assert_eq!(&fasta.name(n)?, name);
            assert_eq!(&fasta.get(n)?.name, name);
        }

        Ok(())
    }

    macro_rules! database_par_iter_test {
        ($name:ident, $func:expr, $threads:expr) => {
            #[test]
            fn $name() -> anyhow::Result<()> {
                rayon::ThreadPoolBuilder::new()
                    .num_threads($threads)
                    .build()
                    .unwrap()
                    .install(|| -> anyhow::Result<()> { $func() })
            }
        };
    }

    fn test_parse_fasta() -> anyhow::Result<()> {
        let path = Path::new(env!("CARGO_MANIFEST_DIR")).join("../fixtures/target.fa");
        let database = Fasta::from_path(&path)?;
        let _ = database.par_iter().collect::<anyhow::Result<Vec<_>>>()?;
        Ok(())
    }

    database_par_iter_test!(test_fasta_par_iter_parse_1_thread, test_parse_fasta, 1);
    database_par_iter_test!(test_fasta_par_iter_parse_2_threads, test_parse_fasta, 2);
    database_par_iter_test!(test_fasta_par_iter_parse_4_threads, test_parse_fasta, 4);
    database_par_iter_test!(test_fasta_par_iter_parse_8_threads, test_parse_fasta, 8);
    database_par_iter_test!(test_fasta_par_iter_parse_16_threads, test_parse_fasta, 16);
    database_par_iter_test!(test_fasta_par_iter_parse_32_threads, test_parse_fasta, 32);
}
