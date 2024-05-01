use std::sync::{Arc, Mutex};

use libnail::structs::{Hmm, Profile, Sequence};
use rayon::iter::{
    plumbing::{bridge, Consumer, Producer, ProducerCallback, UnindexedConsumer},
    IndexedParallelIterator, IntoParallelIterator, ParallelIterator,
};

pub trait Database<I, O> {
    fn iter(&self) -> DbIter<I, O>;
    fn par_iter(&self) -> DbParIter<I, O>;
}

pub struct ProfileCollection {
    profiles: Vec<Arc<Mutex<Profile>>>,
}

impl ProfileCollection {
    pub fn new(profiles: Vec<Profile>) -> Self {
        Self {
            profiles: profiles
                .into_iter()
                .map(|p| Arc::new(Mutex::new(p)))
                .collect(),
        }
    }
}

pub struct SequenceCollection {
    sequences: Vec<Arc<Sequence>>,
}

impl SequenceCollection {
    pub fn new(sequences: Vec<Sequence>) -> Self {
        Self {
            sequences: sequences.into_iter().map(Arc::new).collect(),
        }
    }
}

impl Database<Arc<Mutex<Profile>>, Arc<Mutex<Profile>>> for ProfileCollection {
    fn iter(&self) -> DbIter<Arc<Mutex<Profile>>, Arc<Mutex<Profile>>> {
        DbIter {
            inner: &self.profiles,
            callback: &|p| p.clone(),
        }
    }

    fn par_iter(&self) -> DbParIter<Arc<Mutex<Profile>>, Arc<Mutex<Profile>>> {
        DbParIter {
            inner: &self.profiles,
            callback: &|p| p.clone(),
        }
    }
}

impl Database<Arc<Sequence>, Arc<Sequence>> for SequenceCollection {
    fn iter(&self) -> DbIter<Arc<Sequence>, Arc<Sequence>> {
        DbIter {
            inner: &self.sequences,
            callback: &|s| s.clone(),
        }
    }

    fn par_iter(&self) -> DbParIter<Arc<Sequence>, Arc<Sequence>> {
        DbParIter {
            inner: &self.sequences,
            callback: &|s| s.clone(),
        }
    }
}

impl Database<Arc<Sequence>, Arc<Mutex<Profile>>> for SequenceCollection {
    fn iter(&self) -> DbIter<Arc<Sequence>, Arc<Mutex<Profile>>> {
        DbIter {
            inner: &self.sequences,
            callback: &|seq| {
                let hmm = Hmm::from_blosum_62_and_sequence(seq).unwrap();
                let mut profile = Profile::new(&hmm);
                profile.calibrate_tau(200, 100, 0.04);
                Arc::new(Mutex::new(profile))
            },
        }
    }

    fn par_iter(&self) -> DbParIter<Arc<Sequence>, Arc<Mutex<Profile>>> {
        DbParIter {
            inner: &self.sequences,
            callback: &|seq| {
                let hmm = Hmm::from_blosum_62_and_sequence(seq).unwrap();
                let mut profile = Profile::new(&hmm);
                profile.calibrate_tau(200, 100, 0.04);
                Arc::new(Mutex::new(profile))
            },
        }
    }
}

///
///
///
pub struct DbIter<'a, I, O> {
    inner: &'a [I],
    callback: &'a (dyn Fn(&I) -> O),
}

///
///
///
struct DbProducer<'a, I, O> {
    inner: &'a [I],
    callback: &'a (dyn Fn(&I) -> O + Send + Sync),
}

///
///
///
pub struct DbParIter<'a, I, O> {
    inner: &'a [I],
    callback: &'a (dyn Fn(&I) -> O + Send + Sync),
}

impl<'a, I, O> Iterator for DbIter<'a, I, O> {
    type Item = O;

    fn next(&mut self) -> Option<Self::Item> {
        match self.inner.first() {
            Some(item) => {
                self.inner = &self.inner[1..];
                Some((self.callback)(item))
            }
            None => None,
        }
    }

    // implementing size_hint to always return the
    // exact size of the iterator is REQUIRED for
    // the ExactSizerIterator trait to work; the
    // default impl of size_hint returns (0, None)
    fn size_hint(&self) -> (usize, Option<usize>) {
        let size = self.inner.len();
        (size, Some(size))
    }
}

// rayon expects the the iterators it
// uses to implement DoubleEndedIterator
impl<'a, I, O> DoubleEndedIterator for DbIter<'a, I, O> {
    fn next_back(&mut self) -> Option<Self::Item> {
        match self.inner.last() {
            Some(item) => {
                self.inner = &self.inner[0..self.len() - 1];
                Some((self.callback)(item))
            }
            None => None,
        }
    }
}

impl<'a, I, O> ExactSizeIterator for DbIter<'a, I, O> {}

// ------------
// rayon traits
// ------------

//
//
//
impl<'a, I, O> Producer for DbProducer<'a, I, O>
where
    I: Send + Sync,
    O: Send + Sync,
{
    type Item = O;

    type IntoIter = DbIter<'a, I, O>;

    fn into_iter(self) -> Self::IntoIter {
        DbIter {
            inner: self.inner,
            callback: self.callback,
        }
    }

    fn split_at(self, index: usize) -> (Self, Self) {
        let (left, right) = self.inner.split_at(index);
        (
            Self {
                inner: left,
                callback: self.callback,
            },
            Self {
                inner: right,
                callback: self.callback,
            },
        )
    }
}

//
//
//
impl<'a, I, O> ParallelIterator for DbParIter<'a, I, O>
where
    I: Send + Sync,
    O: Send + Sync,
{
    type Item = O;

    fn drive_unindexed<C>(self, consumer: C) -> C::Result
    where
        C: UnindexedConsumer<Self::Item>,
    {
        bridge(self, consumer)
    }
}

//
//
//
impl<'a, I, O> IndexedParallelIterator for DbParIter<'a, I, O>
where
    I: Send + Sync,
    O: Send + Sync,
{
    fn len(&self) -> usize {
        self.inner.len()
    }

    fn drive<C: Consumer<Self::Item>>(self, consumer: C) -> C::Result {
        bridge(self, consumer)
    }

    fn with_producer<CB: ProducerCallback<Self::Item>>(self, callback: CB) -> CB::Output {
        let producer = DbProducer {
            inner: self.inner,
            callback: self.callback,
        };

        callback.callback(producer)
    }
}

//
//
//
impl<'a, I, O> IntoParallelIterator for &'a dyn Database<I, O>
where
    I: Send + Sync,
    O: Send + Sync,
{
    type Iter = DbParIter<'a, I, O>;

    type Item = O;

    fn into_par_iter(self) -> Self::Iter {
        self.par_iter()
    }
}
