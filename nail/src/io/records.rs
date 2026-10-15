use std::{
    borrow::Cow,
    collections::{hash_map::DefaultHasher, HashSet},
    hash::{Hash, Hasher},
    path::Path,
};

use anyhow::{bail, Context};
use libsail::{
    format::Format,
    index::{Framed, Index},
};
use rayon::prelude::*;

/// An input file with every record framed and named.
pub struct Records(Framed);

impl Records {
    pub fn open(path: impl AsRef<Path>, format: Format) -> anyhow::Result<Self> {
        let path = path.as_ref();

        let mut framed = Framed::open(path, format)
            .with_context(|| format!("failed to open {}", path.display()))?;
        framed.name_index()?;

        // **note: do not remove this check. nail tells
        //         records apart by position, so a repeated
        //         name would run without error. a reader of
        //         the tabular output could not tell two such
        //         records apart, and nothing in the output
        //         would show the problem
        refuse_repeated_names(framed.index(), path)?;

        Ok(Self(framed))
    }

    pub fn len(&self) -> usize {
        self.0.len()
    }

    pub fn name(&self, n: usize) -> anyhow::Result<String> {
        match self.0.index().name(n)? {
            Some(name) => Ok(String::from_utf8(name.into_owned())?),
            None => bail!("record {n} has no name"),
        }
    }

    /// The position of the record named `name`.
    pub fn position(&self, name: &str) -> anyhow::Result<Option<usize>> {
        Ok(self.0.index().position(name.as_bytes())?)
    }

    /// Parse the `n`th record's bytes with `parse`.
    pub fn parse<T>(&self, n: usize, parse: fn(&[u8]) -> anyhow::Result<T>) -> anyhow::Result<T> {
        let bytes = self
            .0
            .record(n)?
            .with_context(|| format!("no record {n}"))?;

        parse(&bytes)
    }
}

/// The `n`th record's name, which must exist and be utf-8.
fn checked_name(index: &Index, n: usize) -> anyhow::Result<Cow<'_, [u8]>> {
    let Some(name) = index.name(n)? else {
        bail!("record {n} has no name");
    };

    // the output prints every name, so a name that is
    // not utf-8 fails here and not after the search
    std::str::from_utf8(&name)?;

    Ok(name)
}

fn refuse_repeated_names(index: &Index, path: &Path) -> anyhow::Result<()> {
    let shards = rayon::current_num_threads();

    let hashes = (0..index.len())
        .into_par_iter()
        .map(|n| {
            let mut hasher = DefaultHasher::new();
            checked_name(index, n)?.hash(&mut hasher);
            Ok(hasher.finish())
        })
        .collect::<anyhow::Result<Vec<u64>>>()?;

    // what: each shard checks the records whose name hash
    //       falls in it, so a repeated name lands twice in
    //       the same shard
    // why: on 2,455,940 targets at 8 threads this takes
    //      0.05 s, against 0.15 s for one HashSet on one
    //      thread and 2.4 s for a binary search per record
    //      over the name index. the gap should grow with
    //      the number of targets and the number of threads,
    //      though no larger run has measured it yet
    (0..shards).into_par_iter().try_for_each(|shard| {
        let mut seen = HashSet::with_capacity(index.len() / shards);

        for (n, &hash) in hashes.iter().enumerate() {
            if hash as usize % shards != shard {
                continue;
            }

            let name = checked_name(index, n)?;

            // the name borrows from the index, so the clone
            // copies a pointer and not the name
            if !seen.insert(name.clone()) {
                bail!(
                    "{} holds more than one record named \"{}\"",
                    path.display(),
                    String::from_utf8_lossy(&name)
                );
            }
        }

        Ok(())
    })
}

#[cfg(test)]
mod tests {
    use std::path::Path;

    use libsail::format::Format;

    use super::*;

    #[test]
    fn a_repeated_name_is_refused() -> anyhow::Result<()> {
        // records 1 and 3 are both named "b"
        let path = Path::new(env!("CARGO_MANIFEST_DIR")).join("../fixtures/repeated-name.fa");

        // one thread puts every name in one shard
        for threads in [1, 2, 8] {
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(threads)
                .build()?;

            let Err(err) = pool.install(|| Records::open(&path, Format::Fasta)) else {
                panic!("a repeated name was accepted at {threads} threads");
            };
            assert!(err
                .to_string()
                .ends_with("holds more than one record named \"b\""));
        }

        Ok(())
    }
}
