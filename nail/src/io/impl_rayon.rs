use crate::io::{
    RecordParser, {Database, DatabaseIter},
};

use rayon::iter::{
    plumbing::{bridge, Consumer, Producer, ProducerCallback, UnindexedConsumer},
    IndexedParallelIterator, IntoParallelIterator, ParallelIterator,
};

impl<R> Database<R>
where
    R: RecordParser,
{
    pub fn par_iter(&'_ self) -> DatabaseParIter<'_, R> {
        DatabaseParIter {
            inner: self.clone(),
            names: self.index.keys().collect(),
        }
    }
}

impl<'a, R> IntoParallelIterator for &'a Database<R>
where
    R: RecordParser,
{
    type Iter = DatabaseParIter<'a, R>;
    type Item = anyhow::Result<R::Record>;

    fn into_par_iter(self) -> Self::Iter {
        self.par_iter()
    }
}

/////

pub struct DatabaseParIter<'a, R>
where
    R: RecordParser,
{
    pub(super) inner: Database<R>,
    pub(super) names: Vec<&'a str>,
}

impl<'a, R> ParallelIterator for DatabaseParIter<'a, R>
where
    R: RecordParser,
{
    type Item = anyhow::Result<R::Record>;

    fn drive_unindexed<C>(self, consumer: C) -> C::Result
    where
        C: UnindexedConsumer<Self::Item>,
    {
        bridge(self, consumer)
    }
}

impl<'a, R> IndexedParallelIterator for DatabaseParIter<'a, R>
where
    R: RecordParser,
{
    fn len(&self) -> usize {
        self.inner.len()
    }

    fn drive<C: Consumer<Self::Item>>(self, consumer: C) -> C::Result {
        bridge(self, consumer)
    }

    fn with_producer<CB: ProducerCallback<Self::Item>>(self, callback: CB) -> CB::Output {
        let producer = DatabaseProducer {
            inner: self.inner,
            names: &self.names,
        };

        callback.callback(producer)
    }
}

pub struct DatabaseProducer<'a, R>
where
    R: RecordParser,
{
    inner: Database<R>,
    names: &'a [&'a str],
}

impl<'a, R> Producer for DatabaseProducer<'a, R>
where
    R: RecordParser,
{
    type Item = anyhow::Result<R::Record>;

    // note: at the time of writing, opaque types are unstable in associated
    //       types e.g. something like this isn't allowed here:
    //         type IntoIter = DatabaseIter<'a, R, impl DoubleEndedIter<Item = &str>>;
    //
    //       but, even if this feature is stablized, it might not fly since:
    //         names.iter() seems to produce Iter<Item = &&str>
    type IntoIter = DatabaseIter<'a, R, std::iter::Copied<std::slice::Iter<'a, &'a str>>>;

    fn into_iter(self) -> Self::IntoIter {
        DatabaseIter {
            inner: self.inner,
            // note: copied() is needed here, since
            // we have &[&str] instead of Vec<&Str>
            keys: self.names.iter().copied(),
        }
    }

    fn split_at(self, index: usize) -> (Self, Self) {
        let (left, right) = self.names.split_at(index);
        (
            DatabaseProducer {
                inner: self.inner.clone(),
                names: left,
            },
            DatabaseProducer {
                inner: self.inner,
                names: right,
            },
        )
    }
}

#[cfg(test)]
mod tests {
    use std::path::Path;

    use crate::io::{Fasta, P7Hmm};

    use rayon::iter::ParallelIterator;

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
