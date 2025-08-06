dyn_clone::clone_trait_object!(<T> Database<T>);
pub trait Database<T>: dyn_clone::DynClone + Send + Sync {
    fn get(&mut self, name: &str) -> Option<T>;
    fn len(&self) -> usize;
    fn iter(&self) -> DatabaseIter<T>;
}

pub struct DatabaseIter<'a, T> {
    pub(crate) inner: Box<dyn Database<T>>,
    pub(crate) names_iter: Box<dyn DoubleEndedIterator<Item = &'a str> + 'a>,
}

impl<'a, T> DatabaseIter<'a, T> {
    pub fn names(self) -> Vec<&'a str> {
        self.names_iter.collect()
    }
}

impl<'a, T> Iterator for DatabaseIter<'a, T> {
    type Item = T;

    fn next(&mut self) -> Option<Self::Item> {
        self.names_iter.next().and_then(|n| self.inner.get(n))
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
impl<'a, T> DoubleEndedIterator for DatabaseIter<'a, T> {
    fn next_back(&mut self) -> Option<Self::Item> {
        self.names_iter.next_back().and_then(|n| self.inner.get(n))
    }
}

impl<'a, T> ExactSizeIterator for DatabaseIter<'a, T> {}
