dyn_clone::clone_trait_object!(<T> Database<T>);
pub trait Database<T>: dyn_clone::DynClone + Send + Sync {
    fn get(&mut self, name: &str) -> Option<T>;
    fn len(&self) -> usize;
    fn iter(&'_ self) -> DatabaseIter<'_, T>;
    fn values(&'_ self) -> DatabaseValues<'_, T>;
}

pub struct DatabaseIter<'a, T> {
    pub(crate) inner: Box<dyn Database<T>>,
    pub(crate) names_iter: Box<dyn DoubleEndedIterator<Item = &'a str> + 'a>,
}

impl<'a, T> Iterator for DatabaseIter<'a, T> {
    type Item = (&'a str, T);

    fn next(&mut self) -> Option<Self::Item> {
        let n = self.names_iter.next()?;
        let v = self.inner.get(n)?;
        Some((n, v))
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

pub struct DatabaseValues<'a, T> {
    pub(crate) inner: Box<dyn Database<T>>,
    pub(crate) names_iter: Box<dyn DoubleEndedIterator<Item = &'a str> + 'a>,
}

impl<'a, T> DatabaseValues<'a, T> {
    pub fn names(self) -> Vec<&'a str> {
        self.names_iter.collect()
    }
}

impl<'a, T> Iterator for DatabaseValues<'a, T> {
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
impl<'a, T> DoubleEndedIterator for DatabaseValues<'a, T> {
    fn next_back(&mut self) -> Option<Self::Item> {
        self.names_iter.next_back().and_then(|n| self.inner.get(n))
    }
}

impl<'a, T> ExactSizeIterator for DatabaseValues<'a, T> {}

#[allow(dead_code)]
mod tmp {
    use std::{
        fs::File,
        io::{BufWriter, Read, Seek, SeekFrom, Write},
        path::{Path, PathBuf},
        sync::{Arc, Mutex},
    };

    use crate::io::{
        Delimiter, DelimiterLocation, Offset, ReadSeekExt, ReadState, RecordParser, SeekableTake,
    };

    use anyhow::{anyhow, bail, Context};
    use indexmap::IndexMap;
    use rayon::iter::{IntoParallelRefMutIterator, ParallelIterator};

    #[derive(Debug, Clone)]
    struct IndexCollisionError {
        key: String,
        first: Offset,
        second: Offset,
    }

    impl std::fmt::Display for IndexCollisionError {
        fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
            write!(
                f,
                "collision of index key: {}\nfirst: {}->{} bytes\nsecond: {}->{} bytes",
                self.key,
                self.first.start,
                self.first.n_bytes,
                self.second.start,
                self.second.n_bytes
            )
        }
    }

    impl std::error::Error for IndexCollisionError {}

    pub struct DatabaseIter<'a, R>
    where
        R: RecordParser,
    {
        pub(crate) inner: Database<R>,
        pub(crate) keys: indexmap::map::Keys<'a, String, Offset>,
    }

    impl<'a, R> Iterator for DatabaseIter<'a, R>
    where
        R: RecordParser,
    {
        type Item = anyhow::Result<R::Record>;

        fn next(&mut self) -> Option<Self::Item> {
            // if the keys are empty, the
            // iterator should terminate
            let key = self.keys.next()?;

            match self.inner.get(key) {
                Some(maybe_rec) => match maybe_rec {
                    // the happy path
                    Ok(rec) => Some(Ok(rec)),
                    // this means parsing the record somehow failed
                    Err(e) => Some(Err(e)),
                },
                // if this happens, the key is somehow not in the
                // index, which shouldn't be possible since
                // the keys iterator is derived from the index map
                None => Some(Err(anyhow!("very strange: index key returned None"))),
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

    pub struct Database<R>
    where
        R: RecordParser,
    {
        pub(crate) path: PathBuf,
        pub(crate) file: File,
        pub(crate) index: Arc<IndexMap<String, Offset>>,
        parser: R,
        buffer: Vec<u8>,
    }

    impl<R> Clone for Database<R>
    where
        R: RecordParser,
    {
        fn clone(&self) -> Self {
            Self {
                path: self.path.clone(),
                file: self.file.try_clone().unwrap(),
                index: self.index.clone(),
                parser: self.parser.clone(),
                buffer: vec![],
            }
        }
    }

    impl<R> Database<R>
    where
        R: RecordParser,
    {
        pub fn from_path<P: AsRef<Path>>(path: P) -> anyhow::Result<Self> {
            let index = Self::build_index::<15, 1>(path.as_ref())?;
            let file = File::open(path.as_ref())?;

            Ok(Self {
                file,
                index: Arc::new(index),
                buffer: Vec::new(),
                parser: R::new(0),
                path: PathBuf::from(path.as_ref()),
            })
        }

        pub fn len(&self) -> usize {
            self.index.len()
        }

        fn iter(&'_ self) -> DatabaseIter<'_, R> {
            DatabaseIter {
                inner: self.clone(),
                keys: self.index.keys(),
            }
        }

        pub fn write<W: Write>(&self, out: W) -> anyhow::Result<()> {
            let mut out = BufWriter::new(out);
            self.iter().try_for_each(|rec| {
                let rec = rec.context("itertor failed to produce a record")?;
                writeln!(out, "{rec}")
                    .context("failed to write record in call to Database::write()")
            })?;

            Ok(())
        }

        pub fn get(&mut self, name: &str) -> Option<anyhow::Result<R::Record>> {
            let offset = self.index.get(name)?;

            self.buffer.resize(offset.n_bytes, 0u8);

            if let Err(e) = self.file.seek(SeekFrom::Start(offset.start)) {
                return Some(Err(e).context("failed to seek to offset in database file"));
            }

            if let Err(e) = self.file.read_exact(&mut self.buffer) {
                return Some(Err(e).context("failed to read from database file"));
            };

            Some(R::parse(&self.buffer))
        }

        fn build_index<const B: usize, const N: usize>(
            path: impl AsRef<Path>,
        ) -> anyhow::Result<IndexMap<String, Offset>> {
            let mut file = File::open(path.as_ref())?;

            let sz = file.metadata().unwrap().len();
            let chunk_sz = sz / N as u64;
            let mut delim_positions: Vec<u64> = vec![];

            // 2 << (B - 1) == 2^B
            let mut buf = vec![0; 2 << (B - 1)];

            use ReadState::*;
            'outer: for mut start in (1..N).map(|i| i as u64 * chunk_sz) {
                file.seek(SeekFrom::Start(start))?;

                'inner: loop {
                    match file.read_with_state(&mut buf)? {
                        Reading(bytes_read) => match R::delimiter_location(&buf[..bytes_read]) {
                            DelimiterLocation::Present(pos) => {
                                delim_positions.push(start + pos);
                                break 'inner;
                            }
                            DelimiterLocation::Ambiguous(n_back) => {
                                file.seek_relative(-(n_back as i64))?;
                                start += bytes_read as u64 - n_back;
                                continue;
                            }
                            DelimiterLocation::Absent => {
                                start += bytes_read as u64;
                            }
                        },
                        Final(bytes_read) => {
                            if let DelimiterLocation::Present(pos) =
                                R::delimiter_location(&buf[..bytes_read])
                            {
                                delim_positions.push(start + pos);
                            }

                            break 'inner;
                        }
                        Done => break 'outer,
                    }
                }
            }

            use Delimiter::*;
            let (block_starts, block_ends) = match R::DELIM_TYPE {
                Initiating => {
                    let mut block_starts = delim_positions;
                    block_starts.dedup();
                    let mut block_ends: Vec<u64> = block_starts.iter().map(|b| b - 1).collect();

                    block_starts.insert(0, 0);
                    block_ends.push(sz - 1);

                    (block_starts, block_ends)
                }
                Terminating => {
                    let mut block_ends = delim_positions
                        .into_iter()
                        .map(|p| {
                            if R::DELIM != b"\n" {
                                // if we have any delimeter other than newline,
                                // we jump past the newline after the delimiter
                                p + R::DELIM_LEN as u64
                            } else {
                                p
                            }
                        })
                        .collect::<Vec<u64>>();

                    block_ends.retain(|e| *e < sz);
                    block_ends.push(sz - 1);
                    block_ends.dedup();

                    let mut block_starts = vec![0];
                    block_starts
                        .extend(block_ends.iter().take(block_ends.len() - 1).map(|b| b + 1));

                    (block_starts, block_ends)
                }
            };

            let mut files_and_starts: Vec<_> = block_starts
                .into_iter()
                .zip(block_ends)
                .map(|(start, end)| {
                    let file = File::open(path.as_ref())?;
                    let limit = end - start + 1;
                    Ok((SeekableTake::new(file, start, limit)?, start))
                    // file.seek(SeekFrom::Start(start))?;
                    // Ok((file.take(limit), start))
                })
                .collect::<anyhow::Result<_>>()?;

            let index = Arc::new(Mutex::new(IndexMap::new()));

            if let Err(e) = files_and_starts
                .par_iter_mut()
                .try_for_each(|(file, start)| {
                    Self::build_chunk::<_, 5_000, 20>(file, index.clone(), Some(*start))
                })
            {
                match e.downcast_ref::<IndexCollisionError>() {
                    Some(e) => {
                        use ReadState::*;
                        let mut buf = vec![0; 2 << (B - 1)];

                        fn search<R>(
                            file: &mut R,
                            buf: &mut [u8],
                            start: u64,
                            target: u64,
                        ) -> anyhow::Result<u64>
                        where
                            R: ReadSeekExt,
                        {
                            let mut pos = start;
                            let mut line_cnt = 0;
                            file.seek(SeekFrom::Start(start))?;
                            while pos <= target {
                                let buf_slice = match file.read_with_state(buf)? {
                                    Reading(_) => buf,
                                    Final(n) => &buf[0..n],
                                    Done => {
                                        break;
                                    }
                                };

                                for &byte in buf_slice {
                                    if byte == b'\n' {
                                        line_cnt += 1;
                                    }

                                    if pos == target {
                                        return Ok(line_cnt);
                                    }

                                    pos += 1;
                                }
                            }

                            Err(anyhow!(
                                "{}\n{}{}\n{}{}\n{}{}\n{}{}",
                                "search failed",
                                "started at byte: ",
                                start,
                                "targeted byte:   ",
                                target,
                                "bytes read:      ",
                                pos,
                                "newline count:   ",
                                line_cnt,
                            ))
                        }

                        let first_line = search(&mut file, &mut buf, 0, e.first.start)
                        .context("failed to recover first record--this probably means there's a bug in nail")
                        .context(e.to_owned())?;

                        let second_line = search(&mut file, &mut buf, e.first.start, e.second.start)
                        .context("failed to recover second record--this probably means there's a bug in nail")
                        .context(e.to_owned())?;

                        use crate::util::term::*;
                        return Err(anyhow!(
                            "{}{RED}\"{}\"{RESET}\n{}{RED}line {}{RESET}\n{}{RED}line {}{RESET}",
                            "the file contains two records that share the identifier ",
                            e.key,
                            "  the first record appears on  ",
                            first_line,
                            "  the second record appears on ",
                            second_line,
                        ));
                    }
                    None => return Err(e),
                }
            }

            Arc::into_inner(index)
                .ok_or(anyhow!("no index"))?
                .into_inner()
                .ok()
                .context("mutex poisoned")
        }

        fn build_chunk<X: Read + Seek, const N: usize, const B: usize>(
            mut data_chunk: X,
            index: Arc<Mutex<IndexMap<String, Offset>>>,
            start: Option<u64>,
        ) -> anyhow::Result<()> {
            let process_fn = |recs: &mut Vec<(String, Offset)>, done: bool| {
                let full = recs.len() == N;

                let guard = if full || done {
                    match index.lock() {
                        Ok(guard) => Some(guard),
                        Err(_) => bail!("mutex poisoned"),
                    }
                } else {
                    use std::sync::TryLockError::*;
                    match index.try_lock() {
                        Ok(guard) => Some(guard),
                        Err(WouldBlock) => None,
                        Err(Poisoned(_)) => bail!("mutex poisoned"),
                    }
                };

                match guard {
                    Some(mut guard) => {
                        guard.reserve(recs.len());
                        for (key, offset) in recs.drain(..) {
                            use indexmap::map::Entry;

                            match guard.entry(key) {
                                Entry::Vacant(e) => {
                                    e.insert(offset);
                                }
                                Entry::Occupied(e) => {
                                    return Err(IndexCollisionError {
                                        key: e.key().clone(),
                                        first: e.get().clone(),
                                        second: offset,
                                    }
                                    .into());
                                }
                            }
                        }

                        Ok(())
                    }

                    None => Ok(()),
                }
            };

            use ReadState::*;
            let mut pos = start.unwrap_or(0);
            let mut parser = R::new(pos);
            let mut on_first_record = true;
            let mut buf = vec![0; 2 << (B - 1)];
            let mut names_and_offsets: Vec<(String, Offset)> = Vec::with_capacity(N);

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
                        if let Delimiter::Initiating = R::DELIM_TYPE {
                            if *buf.last().unwrap() != b'\n' {
                                buf.push(b'\n');
                            }
                            buf.extend(R::DELIM);
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
                        if let (true, Delimiter::Initiating) = (on_first_record, R::DELIM_TYPE) {
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
}

mod index {}
