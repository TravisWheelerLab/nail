use std::path::Path;

use anyhow::{bail, Context};
use libsail::{format::Format, index::Framed, source::MmapSource};

/// An input file mapped into memory, with every record framed and named.
pub struct Records(Framed<MmapSource>);

impl Records {
    pub fn map(path: impl AsRef<Path>, format: Format) -> anyhow::Result<Self> {
        let path = path.as_ref();

        // SAFETY: nothing may modify or truncate the file
        //         while it is mapped. nail never writes to
        //         its inputs
        let mut framed = unsafe { Framed::map(path, format) }
            .with_context(|| format!("failed to map {}", path.display()))?;
        framed.name_index()?;

        let records = Self(framed);

        // the index this replaced refused a repeated name,
        // and a lookup by name returns only the first record
        for n in 0..records.len() {
            let name = records.name(n)?;
            let first = records.position(&name)?;
            if first != Some(n) {
                bail!(
                    "{} holds more than one record named \"{name}\"",
                    path.display()
                );
            }
        }

        Ok(records)
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
