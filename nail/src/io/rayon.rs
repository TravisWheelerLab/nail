use libnail::structs::Sequence;
use rayon::iter::{
    plumbing::{bridge, Consumer, Producer, ProducerCallback, UnindexedConsumer},
    IndexedParallelIterator, IntoParallelIterator, ParallelIterator,
};

use super::{Fasta, SequenceDatabase, SequenceDatabaseIter};

pub struct SequenceDatabaseParIter<'a> {
    pub(super) inner: Box<dyn SequenceDatabase>,
    pub(super) names: Vec<&'a str>,
}

impl<'a> ParallelIterator for SequenceDatabaseParIter<'a> {
    type Item = Sequence;

    fn drive_unindexed<C>(self, consumer: C) -> C::Result
    where
        C: UnindexedConsumer<Self::Item>,
    {
        bridge(self, consumer)
    }
}

impl<'a> IndexedParallelIterator for SequenceDatabaseParIter<'a> {
    fn len(&self) -> usize {
        self.inner.len()
    }

    fn drive<C: Consumer<Self::Item>>(self, consumer: C) -> C::Result {
        bridge(self, consumer)
    }

    fn with_producer<CB: ProducerCallback<Self::Item>>(self, callback: CB) -> CB::Output {
        let producer = SequenceDatabaseProducer {
            inner: self.inner,
            names: &self.names,
        };

        callback.callback(producer)
    }
}

impl Fasta {
    fn par_iter(&self) -> SequenceDatabaseParIter {
        SequenceDatabaseParIter {
            inner: Box::new(self.clone()),
            names: self.index.offsets.keys().map(|s| s.as_str()).collect(),
        }
    }
}

impl<'a> IntoParallelIterator for &'a Fasta {
    type Iter = SequenceDatabaseParIter<'a>;

    type Item = Sequence;

    fn into_par_iter(self) -> Self::Iter {
        self.par_iter()
    }
}

pub struct SequenceDatabaseProducer<'a> {
    inner: Box<dyn SequenceDatabase>,
    names: &'a [&'a str],
}

impl<'a> Producer for SequenceDatabaseProducer<'a> {
    type Item = Sequence;

    type IntoIter = SequenceDatabaseIter<'a>;

    fn into_iter(self) -> Self::IntoIter {
        SequenceDatabaseIter {
            inner: self.inner,
            names_iter: Box::new(self.names.iter().copied()),
        }
    }

    fn split_at(self, index: usize) -> (Self, Self) {
        let (left, right) = self.names.split_at(index);
        (
            SequenceDatabaseProducer {
                inner: self.inner.clone(),
                names: left,
            },
            SequenceDatabaseProducer {
                inner: self.inner,
                names: right,
            },
        )
    }
}
