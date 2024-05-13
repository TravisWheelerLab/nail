use std::{
    collections::HashMap,
    sync::{Arc, Mutex},
};

use libnail::structs::{Hmm, Profile, Sequence};
use rayon::iter::{
    plumbing::{bridge, Consumer, Producer, ProducerCallback, UnindexedConsumer},
    IndexedParallelIterator, IntoParallelIterator, ParallelIterator,
};

use crate::id::{next_id, Id};

pub trait Database<O> {
    fn iter(&self) -> DbIter<O>;
    fn par_iter(&self) -> DbParIter<O>;
    fn get(&self, id: &Id) -> Option<O>;
}

pub struct ProfileCollection {
    profile_map: HashMap<Id, Arc<Mutex<Profile>>>,
    ids: Vec<Id>,
}

impl ProfileCollection {
    pub fn new(profiles: Vec<Profile>) -> Self {
        let ids: Vec<_> = profiles.iter().map(|_| next_id()).collect();
        let mut profile_map: HashMap<_, _> = HashMap::new();

        profiles
            .into_iter()
            .map(|p| Arc::new(Mutex::new(p)))
            .zip(&ids)
            .for_each(|(p, &i)| {
                profile_map.insert(i, p);
            });

        Self { profile_map, ids }
    }

    pub fn len(&self) -> usize {
        self.profile_map.len()
    }

    pub fn ids(&self) -> &[Id] {
        &self.ids
    }
}

pub struct SequenceCollection {
    sequence_map: HashMap<Id, Arc<Sequence>>,
    ids: Vec<Id>,
}

impl SequenceCollection {
    pub fn new(sequences: Vec<Sequence>) -> Self {
        let ids: Vec<_> = sequences.iter().map(|_| next_id()).collect();
        let mut sequence_map: HashMap<_, _> = HashMap::new();

        sequences
            .into_iter()
            .map(Arc::new)
            .zip(&ids)
            .for_each(|(p, &i)| {
                sequence_map.insert(i, p);
            });
        Self { sequence_map, ids }
    }

    pub fn len(&self) -> usize {
        self.sequence_map.len()
    }

    pub fn ids(&self) -> &[Id] {
        &self.ids
    }
}

impl Database<Arc<Mutex<Profile>>> for ProfileCollection {
    fn iter(&self) -> DbIter<Arc<Mutex<Profile>>> {
        DbIter {
            inner: self,
            ids: &self.ids,
        }
    }

    fn par_iter(&self) -> DbParIter<Arc<Mutex<Profile>>> {
        DbParIter {
            inner: self,
            ids: &self.ids,
        }
    }

    fn get(&self, id: &Id) -> Option<Arc<Mutex<Profile>>> {
        self.profile_map.get(id).cloned()
    }
}

impl Database<Arc<Sequence>> for SequenceCollection {
    fn iter(&self) -> DbIter<Arc<Sequence>> {
        DbIter {
            inner: self,
            ids: &self.ids,
        }
    }

    fn par_iter(&self) -> DbParIter<Arc<Sequence>> {
        DbParIter {
            inner: self,
            ids: &self.ids,
        }
    }

    fn get(&self, id: &Id) -> Option<Arc<Sequence>> {
        self.sequence_map.get(id).cloned()
    }
}

impl Database<Arc<Mutex<Profile>>> for SequenceCollection {
    fn iter(&self) -> DbIter<Arc<Mutex<Profile>>> {
        DbIter {
            inner: self,
            ids: &self.ids,
        }
    }

    fn par_iter(&self) -> DbParIter<Arc<Mutex<Profile>>> {
        DbParIter {
            inner: self,
            ids: &self.ids,
        }
    }

    fn get(&self, id: &Id) -> Option<Arc<Mutex<Profile>>> {
        match self.sequence_map.get(id) {
            Some(seq) => {
                let hmm = Hmm::from_blosum_62_and_sequence(seq)
                    .expect("failed to create Hmm from Sequence");
                let mut profile = Profile::new(&hmm);
                profile.calibrate_tau(200, 100, 0.04);
                Some(Arc::new(Mutex::new(profile)))
            }
            None => None,
        }
    }
}

pub struct DbIter<'a, O> {
    inner: &'a dyn Database<O>,
    ids: &'a [Id],
}

pub struct DbProducer<'a, O> {
    inner: &'a (dyn Database<O> + Send + Sync),
    ids: &'a [Id],
}

pub struct DbParIter<'a, O> {
    inner: &'a (dyn Database<O> + Send + Sync),
    ids: &'a [Id],
}

impl<'a, O> Iterator for DbIter<'a, O> {
    type Item = O;

    fn next(&mut self) -> Option<Self::Item> {
        match self.ids.first() {
            Some(id) => {
                self.ids = &self.ids[1..];
                match self.inner.get(id) {
                    Some(item) => Some(item),
                    // being very cautious here
                    None => panic!("invalid Id in DbIter"),
                }
            }
            None => None,
        }
    }

    // implementing size_hint to always return the
    // exact size of the iterator is REQUIRED for
    // the ExactSizerIterator trait to work; the
    // default impl of size_hint returns (0, None)
    fn size_hint(&self) -> (usize, Option<usize>) {
        let size = self.ids.len();
        (size, Some(size))
    }
}

// rayon expects the the iterators it
// uses to implement DoubleEndedIterator
impl<'a, O> DoubleEndedIterator for DbIter<'a, O> {
    fn next_back(&mut self) -> Option<Self::Item> {
        match self.ids.last() {
            Some(id) => {
                self.ids = &self.ids[0..self.len() - 1];
                match self.inner.get(id) {
                    Some(item) => Some(item),
                    None => panic!("invalid Id in DbIter"),
                }
            }
            None => None,
        }
    }
}

impl<'a, O> ExactSizeIterator for DbIter<'a, O> {}

// ------------
// rayon traits
// ------------

//
//
//

impl<'a, O> Producer for DbProducer<'a, O>
where
    O: Send + Sync,
{
    type Item = O;

    type IntoIter = DbIter<'a, O>;

    fn into_iter(self) -> Self::IntoIter {
        DbIter {
            inner: self.inner,
            ids: self.ids,
        }
    }

    fn split_at(self, index: usize) -> (Self, Self) {
        let (left, right) = self.ids.split_at(index);
        (
            Self {
                inner: self.inner,
                ids: left,
            },
            Self {
                inner: self.inner,
                ids: right,
            },
        )
    }
}

impl<'a, O> ParallelIterator for DbParIter<'a, O>
where
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

impl<'a, O> IndexedParallelIterator for DbParIter<'a, O>
where
    O: Send + Sync,
{
    fn len(&self) -> usize {
        self.ids.len()
    }

    fn drive<C: Consumer<Self::Item>>(self, consumer: C) -> C::Result {
        bridge(self, consumer)
    }

    fn with_producer<CB: ProducerCallback<Self::Item>>(self, callback: CB) -> CB::Output {
        let producer = DbProducer {
            inner: self.inner,
            ids: self.ids,
        };

        callback.callback(producer)
    }
}

impl<'a, O> IntoParallelIterator for &'a dyn Database<O>
where
    O: Send + Sync,
{
    type Iter = DbParIter<'a, O>;

    type Item = O;

    fn into_par_iter(self) -> Self::Iter {
        self.par_iter()
    }
}
