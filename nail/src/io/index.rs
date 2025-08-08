use std::{
    fs::File,
    io::{Read, Seek, SeekFrom},
    marker::PhantomData,
    path::Path,
    sync::{Arc, Mutex},
};

use anyhow::{anyhow, bail, Context};
use indexmap::IndexMap;
use rayon::iter::{IntoParallelRefMutIterator, ParallelIterator};

use crate::io::{ReadSeekExt, SeekableTake};

use super::ReadState;

pub enum Delimiter {
    Initiating,
    Terminating,
}

pub trait RecordParser: Send + Sync + 'static {
    const DELIM: &'static [u8];
    const DELIM_TYPE: Delimiter;
    type Offset;
    type Record;

    fn new() -> Self;
    fn offset(&mut self, line: &[u8], start: u64) -> Option<(String, Self::Offset)>;
    fn parse(buf: &[u8]) -> anyhow::Result<Self::Record>;
}

pub trait IndexInner<O>: Send + Sync + 'static {
    fn new() -> Self;
    fn len(&self) -> usize;
    fn get(&self, name: &str) -> Option<&O>;
    fn extend<T>(&mut self, iter: T)
    where
        T: IntoIterator<Item = (String, O)>;
}

impl<O> IndexInner<O> for IndexMap<String, O>
where
    O: Send + Sync + 'static,
{
    fn new() -> Self {
        Self::new()
    }

    fn len(&self) -> usize {
        self.len()
    }

    fn get(&self, name: &str) -> Option<&O> {
        self.get(name)
    }

    fn extend<T>(&mut self, iter: T)
    where
        T: IntoIterator<Item = (String, O)>,
    {
        // automatic delegation fails here since
        // Extend is actually a trait in std::iter
        std::iter::Extend::extend(self, iter);
    }
}

pub struct Index<I, P>
where
    I: IndexInner<P::Offset>,
    P: RecordParser,
{
    pub(crate) inner: I,
    _marker: std::marker::PhantomData<P>,
}

impl<I, P> Index<I, P>
where
    P: RecordParser,
    I: IndexInner<P::Offset>,
{
    fn new() -> Self {
        Self {
            inner: I::new(),
            _marker: PhantomData,
        }
    }

    pub fn len(&self) -> usize {
        self.inner.len()
    }

    pub fn get(&self, name: &str) -> Option<&P::Offset> {
        self.inner.get(name)
    }

    fn extend<T>(&mut self, iter: T)
    where
        T: IntoIterator<Item = (String, P::Offset)>,
    {
        self.inner.extend(iter);
    }

    pub fn from_path<const B: usize, const N: usize>(
        path: impl AsRef<Path>,
    ) -> anyhow::Result<Self> {
        let mut file = File::open(path.as_ref())?;
        let sz = file.metadata().unwrap().len();
        let chunk_sz = sz / N as u64;

        // 2 << (B - 1) == 2^B
        let mut buf = vec![0; 2 << (B - 1)];
        let mut block_starts: Vec<u64> = vec![];

        let (delim_first, delim_rem) = P::DELIM.split_first().ok_or(anyhow!("empty delimiter"))?;

        let delim_len = delim_rem.len();

        use ReadState::*;
        use SeekFrom::{Current, Start};
        'outer: for mut start in (0..N).map(|i| i as u64 * chunk_sz) {
            file.seek(Start(start))?;

            'inner: loop {
                match file.read_with_state(&mut buf)? {
                    Reading(n) => {
                        for (pos, b) in buf[..n].iter().enumerate() {
                            if b == delim_first {
                                let st = pos + 1;
                                let buf_len = buf[st..].len();

                                if buf_len < delim_len {
                                    file.seek(Current(-((buf_len + 1) as i64)))?;
                                    continue;
                                };

                                if &buf[st..st + delim_len] == delim_rem {
                                    block_starts.push(start + pos as u64);
                                    break 'inner;
                                }
                            }
                        }
                        start += n as u64;
                    }
                    Final(n) => {
                        for (pos, b) in buf[..n].iter().enumerate() {
                            if b == delim_first {
                                let st = pos + 1;
                                let buf_len = buf[st..].len();

                                if buf_len < delim_len {
                                    break 'outer;
                                };

                                if &buf[st..st + delim_len] == delim_rem {
                                    block_starts.push(start + st as u64);
                                    break 'inner;
                                }
                            }
                        }
                    }
                    Done => break 'outer,
                }
            }
        }

        let mut block_ends: Vec<u64> = block_starts.iter().skip(1).map(|b| b - 1).collect();
        block_ends.push(sz - 1);

        let mut files_and_starts: Vec<_> = block_starts
            .into_iter()
            .zip(block_ends)
            .map(|(start, end)| {
                let file = File::open(path.as_ref())?;
                let limit = end - start + 1;
                Ok((SeekableTake::new(file, start, limit)?, start))
            })
            .collect::<anyhow::Result<_>>()?;

        let index = Arc::new(Mutex::new(Self::new()));

        files_and_starts
            .par_iter_mut()
            .try_for_each(|(file, start)| {
                Self::build_chunk::<_, 5_000, 20>(file, index.clone(), Some(*start))
            })?;

        Arc::into_inner(index)
            .ok_or(anyhow!("no index"))?
            .into_inner()
            .ok()
            .context("mutex poisoned")
    }

    fn build_chunk<R: Read + Seek, const N: usize, const B: usize>(
        mut data_chunk: R,
        index: Arc<Mutex<Self>>,
        start: Option<u64>,
    ) -> anyhow::Result<()> {
        let process_fn = |recs: &mut Vec<(String, P::Offset)>, done: bool| {
            match index.try_lock() {
                Ok(mut guard) => {
                    guard.extend(recs.drain(..));
                }
                Err(std::sync::TryLockError::WouldBlock) => {
                    if recs.len() == N || done {
                        index
                            .lock()
                            .map_err(|_| anyhow!("mutex poisoned"))?
                            .extend(recs.drain(..));
                    }
                }
                Err(_) => bail!("mutex poisoned"),
            }

            Ok(())
        };

        use ReadState::*;
        let mut pos = start.unwrap_or(0);
        let mut parser = P::new();
        let mut on_first_record = true;
        let mut buf = vec![0; 2 << (B - 1)];
        let mut names_and_offsets: Vec<(String, P::Offset)> = Vec::with_capacity(N);

        while let Ok(read_state) = data_chunk.read_with_state(&mut buf) {
            let buf_slice = match read_state {
                Reading(n) => {
                    let last_newline_pos = buf[..n].iter().rposition(|&b| b == b'\n').ok_or(
                        anyhow!("failed to find any newlines in buf slice of length: {}", n),
                    )?;

                    let n_bytes_back = (n - last_newline_pos - 1) as i64;
                    data_chunk.seek_relative(-n_bytes_back)?;

                    &buf[..=last_newline_pos]
                }
                Final(n) => {
                    buf.resize(n, 0);
                    if let Delimiter::Initiating = P::DELIM_TYPE {
                        if *buf.last().unwrap() != b'\n' {
                            buf.push(b'\n');
                        }
                        buf.extend(P::DELIM);
                    }
                    &buf
                }
                Done => {
                    break;
                }
            };

            for line in buf_slice.split_inclusive(|b| *b == b'\n') {
                if let Some(rec) = parser.offset(line, pos) {
                    names_and_offsets.push(rec);

                    // if this is the first record AND the delimiter preceeds the records
                    // in the file, then the parser will produce an empty record
                    if let (true, Delimiter::Initiating) = (on_first_record, P::DELIM_TYPE) {
                        names_and_offsets.pop();
                    }

                    on_first_record = false;
                    process_fn(&mut names_and_offsets, false)?;
                }

                pos += line.len() as u64;
            }
        }

        process_fn(&mut names_and_offsets, true)?;

        Ok(())
    }
}
