use super::{Database, DatabaseIter, Fasta, P7Hmm};
use libnail::structs::{Profile, Sequence};

use rayon::iter::{
    plumbing::{bridge, Consumer, Producer, ProducerCallback, UnindexedConsumer},
    IndexedParallelIterator, IntoParallelIterator, ParallelIterator,
};

impl Fasta {
    pub fn par_iter(&self) -> DatabaseParIter<Sequence> {
        DatabaseParIter {
            inner: Box::new(self.clone()),
            names: self.index.inner.keys().map(|s| s.as_str()).collect(),
        }
    }
}

impl<'a> IntoParallelIterator for &'a Fasta {
    type Iter = DatabaseParIter<'a, Sequence>;

    type Item = Sequence;

    fn into_par_iter(self) -> Self::Iter {
        self.par_iter()
    }
}

impl P7Hmm {
    pub fn par_iter(&self) -> DatabaseParIter<Profile> {
        DatabaseParIter {
            inner: Box::new(self.clone()),
            names: self.index.offsets.keys().map(|s| s.as_str()).collect(),
        }
    }
}

impl<'a> IntoParallelIterator for &'a P7Hmm {
    type Iter = DatabaseParIter<'a, Profile>;

    type Item = Profile;

    fn into_par_iter(self) -> Self::Iter {
        self.par_iter()
    }
}

/////

pub struct DatabaseParIter<'a, T> {
    pub(super) inner: Box<dyn Database<T>>,
    pub(super) names: Vec<&'a str>,
}

impl<'a, T> ParallelIterator for DatabaseParIter<'a, T>
where
    T: Send + Sync,
{
    type Item = T;

    fn drive_unindexed<C>(self, consumer: C) -> C::Result
    where
        C: UnindexedConsumer<Self::Item>,
    {
        bridge(self, consumer)
    }
}

impl<'a, T> IndexedParallelIterator for DatabaseParIter<'a, T>
where
    T: Send + Sync,
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

pub struct DatabaseProducer<'a, T> {
    inner: Box<dyn Database<T>>,
    names: &'a [&'a str],
}

impl<'a, T> Producer for DatabaseProducer<'a, T> {
    type Item = T;

    type IntoIter = DatabaseIter<'a, T>;

    fn into_iter(self) -> Self::IntoIter {
        DatabaseIter {
            inner: self.inner,
            names_iter: Box::new(self.names.iter().copied()),
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
