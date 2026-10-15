use std::path::Path;

use libnail::structs::Profile;
use libsail::format::Format;
use rayon::prelude::*;

use crate::io::Records;

/// A p7hmm file and where each of its models sits in it.
pub struct P7Hmm(Records);

impl P7Hmm {
    pub fn from_path(path: impl AsRef<Path>) -> anyhow::Result<Self> {
        Ok(Self(Records::open(path, Format::Hmm)?))
    }

    pub fn len(&self) -> usize {
        self.0.len()
    }

    /// The name of the `n`th model.
    pub fn name(&self, n: usize) -> anyhow::Result<String> {
        self.0.name(n)
    }

    /// The `n`th model.
    pub fn get(&self, n: usize) -> anyhow::Result<Profile> {
        self.0.parse(n, Profile::from_p7hmm)
    }

    pub fn iter(&self) -> impl ExactSizeIterator<Item = anyhow::Result<Profile>> + '_ {
        (0..self.len()).map(|n| self.get(n))
    }

    pub fn par_iter(&self) -> impl IndexedParallelIterator<Item = anyhow::Result<Profile>> + '_ {
        (0..self.len()).into_par_iter().map(|n| self.get(n))
    }
}

#[cfg(test)]
mod tests {
    use std::{fs::read_to_string, path::Path};

    use rayon::iter::ParallelIterator;

    use super::*;

    #[test]
    fn test_p7hmm_iter_parse() -> anyhow::Result<()> {
        let path = Path::new(env!("CARGO_MANIFEST_DIR")).join("../fixtures/query.hmm");
        let database = P7Hmm::from_path(&path)?;
        let _ = database.iter().collect::<anyhow::Result<Vec<_>>>()?;
        Ok(())
    }

    #[test]
    fn test_p7hmm_records_in_file_order() -> anyhow::Result<()> {
        let path = Path::new(env!("CARGO_MANIFEST_DIR")).join("../fixtures/query.hmm");

        // the name of each model, by a plain scan of the NAME lines
        let names: Vec<String> = read_to_string(&path)?
            .lines()
            .filter(|line| line.starts_with("NAME"))
            .map(|line| line.split_whitespace().nth(1).unwrap().to_string())
            .collect();

        let p7hmm = P7Hmm::from_path(&path)?;
        assert_eq!(p7hmm.len(), names.len());

        for (n, name) in names.iter().enumerate() {
            assert_eq!(&p7hmm.name(n)?, name);
            assert_eq!(&p7hmm.get(n)?.name, name);
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

    fn test_parse_p7hmm() -> anyhow::Result<()> {
        let path = Path::new(env!("CARGO_MANIFEST_DIR")).join("../fixtures/query.hmm");
        let database = P7Hmm::from_path(&path)?;
        let _ = database.par_iter().collect::<anyhow::Result<Vec<_>>>()?;
        Ok(())
    }

    database_par_iter_test!(test_p7hmm_par_iter_parse_1_thread, test_parse_p7hmm, 1);
    database_par_iter_test!(test_p7hmm_par_iter_parse_2_threads, test_parse_p7hmm, 2);
    database_par_iter_test!(test_p7hmm_par_iter_parse_4_threads, test_parse_p7hmm, 4);
    database_par_iter_test!(test_p7hmm_par_iter_parse_8_threads, test_parse_p7hmm, 8);
    database_par_iter_test!(test_p7hmm_par_iter_parse_16_threads, test_parse_p7hmm, 16);
    database_par_iter_test!(test_p7hmm_par_iter_parse_32_threads, test_parse_p7hmm, 32);
}
