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
            names: self.index.keys().map(|s| s.as_str()).collect(),
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
