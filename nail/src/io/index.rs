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

#[derive(Debug)]
pub enum DelimiterLocation {
    Present(u64),
    Ambiguous(u64),
    Absent,
}

const fn tail(bytes: &'static [u8]) -> &'static [u8] {
    // ASSUMPTION:
    //   bytes.len() >= 1
    unsafe { core::slice::from_raw_parts(bytes.as_ptr().add(1), bytes.len() - 1) }
}

pub trait RecordParser: Send + Sync + 'static {
    const DELIM: &'static [u8];
    const DELIM_HEAD: u8 = Self::DELIM[0];
    const DELIM_TAIL: &'static [u8] = tail(Self::DELIM);
    const DELIM_LEN: usize = Self::DELIM.len();
    const DELIM_TAIL_LEN: usize = Self::DELIM_TAIL.len();
    const DELIM_TYPE: Delimiter;
    type Offset: Offset;
    type Record;

    fn new(start_pos: u64) -> Self;
    fn offset(&mut self, line: &[u8], start: u64) -> Option<(String, Self::Offset)>;
    fn parse(buf: &[u8]) -> anyhow::Result<Self::Record>;

    fn delimiter_location(buf: &[u8]) -> DelimiterLocation {
        match buf.windows(Self::DELIM_LEN).position(|w| w == Self::DELIM) {
            Some(pos) => DelimiterLocation::Present(pos as u64),
            None => {
                let n_skip = buf.len().saturating_sub(Self::DELIM_TAIL_LEN);
                let tail_len = buf.len() - n_skip;
                match buf.iter().skip(n_skip).position(|b| *b == Self::DELIM_HEAD) {
                    Some(pos) => DelimiterLocation::Ambiguous((tail_len - pos - 1) as u64),
                    None => DelimiterLocation::Absent,
                }
            }
        }
    }
}

pub trait Offset: Send + Sync + Clone + 'static {}

pub trait IndexInner<O>: Send + Sync {
    fn new() -> Self;
    fn len(&self) -> usize;
    fn get(&self, name: &str) -> Option<&O>;
    fn extend<T>(&mut self, iter: T) -> anyhow::Result<()>
    where
        T: IntoIterator<Item = (String, O)>;
}

struct IndexCollisionError<O>
where
    O: Offset,
{
    key: String,
    first: O,
    second: O,
}

impl<O> std::fmt::Debug for IndexCollisionError<O>
where
    O: Offset,
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        todo!()
    }
}

impl<O> std::fmt::Display for IndexCollisionError<O>
where
    O: Offset,
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        todo!()
    }
}

impl<O> std::error::Error for IndexCollisionError<O> where O: Offset {}

impl<O> IndexInner<O> for IndexMap<String, O>
where
    O: Offset,
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

    fn extend<T>(&mut self, iter: T) -> anyhow::Result<()>
    where
        T: IntoIterator<Item = (String, O)>,
    {
        let iter = iter.into_iter();

        // note: this resize logic is copied and
        //       modifed from hashbrown in std
        let reserve = if self.is_empty() {
            iter.size_hint().0
        } else {
            iter.size_hint().0.div_ceil(2)
        };

        self.reserve(reserve);

        for (key, offset) in iter {
            use indexmap::map::Entry;

            match self.entry(key) {
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

    fn extend<T>(&mut self, iter: T) -> anyhow::Result<()>
    where
        T: IntoIterator<Item = (String, P::Offset)>,
    {
        self.inner.extend(iter)
    }

    pub fn from_path<const B: usize, const N: usize>(
        path: impl AsRef<Path>,
    ) -> anyhow::Result<Self> {
        let mut file = File::open(path.as_ref())?;

        let sz = file.metadata().unwrap().len();
        let chunk_sz = sz / N as u64;
        let mut delim_positions: Vec<u64> = vec![];

        // 2 << (B - 1) == 2^B
        let mut buf = vec![0; 2 << (B - 1)];

        'outer: for mut start in (1..N).map(|i| i as u64 * chunk_sz) {
            file.seek(SeekFrom::Start(start))?;

            'inner: loop {
                match file.read_with_state(&mut buf)? {
                    ReadState::Reading(bytes_read) => {
                        match P::delimiter_location(&buf[..bytes_read]) {
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
                        }
                    }
                    ReadState::Final(bytes_read) => {
                        if let DelimiterLocation::Present(pos) =
                            P::delimiter_location(&buf[..bytes_read])
                        {
                            delim_positions.push(start + pos);
                        }

                        break 'inner;
                    }
                    ReadState::Done => break 'outer,
                }
            }
        }

        let (block_starts, block_ends) = match P::DELIM_TYPE {
            Delimiter::Initiating => {
                let mut block_starts = delim_positions;
                block_starts.dedup();
                let mut block_ends: Vec<u64> = block_starts.iter().map(|b| b - 1).collect();

                block_starts.insert(0, 0);
                block_ends.push(sz - 1);

                (block_starts, block_ends)
            }
            Delimiter::Terminating => {
                let mut block_ends = delim_positions
                    .into_iter()
                    .map(|p| {
                        if P::DELIM != b"\n" {
                            // if we have any delimeter other than newline,
                            // we jump past the newline after the delimiter
                            p + P::DELIM_LEN as u64
                        } else {
                            p
                        }
                    })
                    .collect::<Vec<u64>>();

                block_ends.retain(|e| *e < sz);
                block_ends.push(sz - 1);
                block_ends.dedup();

                let mut block_starts = vec![0];
                block_starts.extend(block_ends.iter().take(block_ends.len() - 1).map(|b| b + 1));

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
            let full = recs.len() == N;

            let guard = if full || done {
                match index.lock() {
                    Ok(guard) => Some(guard),
                    Err(_) => bail!("mutex poisoned"),
                }
            } else {
                match index.try_lock() {
                    Ok(guard) => Some(guard),
                    Err(std::sync::TryLockError::WouldBlock) => None,
                    Err(std::sync::TryLockError::Poisoned(_)) => bail!("mutex poisoned"),
                }
            };

            match guard {
                Some(mut guard) => match guard.extend(recs.drain(..)) {
                    Ok(_) => Ok(()),
                    Err(e) => match e.downcast::<IndexCollisionError<P::Offset>>() {
                        Ok(e) => {
                            // TODO: figure out a way to dispatch the
                            //       error into a routine that will print
                            //       the line numbers of the collisions
                            Ok(())
                        }
                        Err(_) => bail!("failed to downcast to concrete IndexCollisionError<O>"),
                    },
                },
                None => Ok(()),
            }
        };

        use ReadState::*;
        let mut pos = start.unwrap_or(0);
        let mut parser = P::new(pos);
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

#[cfg(test)]
mod tests {
    use super::*;

    #[derive(Clone)]
    struct DummyOffset {}
    impl Offset for DummyOffset {}

    macro_rules! dummy_parser {
        ($name:ident, $delim:expr) => {
            #[derive(Debug)]
            struct $name;

            impl RecordParser for $name {
                const DELIM: &'static [u8] = $delim;
                const DELIM_TYPE: Delimiter = Delimiter::Initiating;
                type Offset = DummyOffset;
                type Record = ();

                fn new(_: u64) -> Self {
                    $name
                }
                fn offset(&mut self, _: &[u8], _: u64) -> Option<(String, Self::Offset)> {
                    None
                }
                fn parse(_: &[u8]) -> anyhow::Result<Self::Record> {
                    Ok(())
                }
            }
        };
    }

    fn present_at<P: RecordParser>(buf: &[u8], at: u64) {
        match P::delimiter_location(buf) {
            DelimiterLocation::Present(pos) => assert_eq!(pos, at),
            x => panic!("expected Present({at}), got {x:?}"),
        }
    }
    fn ambiguous_at<P: RecordParser>(buf: &[u8], at: u64) {
        match P::delimiter_location(buf) {
            DelimiterLocation::Ambiguous(pos) => assert_eq!(pos, at),
            x => panic!("expected Ambiguous({at}), got {x:?}"),
        }
    }
    fn absent<P: RecordParser>(buf: &[u8]) {
        match P::delimiter_location(buf) {
            DelimiterLocation::Absent => {}
            x => panic!("expected Absent, got {x:?}"),
        }
    }

    // ---- 1-byte delim: ">" ----
    dummy_parser!(One, b">");

    #[test]
    fn one_present() {
        present_at::<One>(b">start", 0);
        present_at::<One>(b"xx>yy", 2);
        present_at::<One>(b"aaa>bbb>ccc", 3);
    }

    #[test]
    fn one_absent_and_no_ambiguous() {
        absent::<One>(b"");
        absent::<One>(b"nope");
        present_at::<One>(b"end>", 3);
    }

    // ---- 2-byte delim: "//" ----
    dummy_parser!(Two, b"//");

    #[test]
    fn two_present() {
        present_at::<Two>(b"//heh", 0);
        present_at::<Two>(b"xx//yy", 2);
        present_at::<Two>(b"/not yet//now", 8);
        present_at::<Two>(b"its/at/the/end//", 14);
    }

    #[test]
    fn two_ambiguous_tail() {
        ambiguous_at::<Two>(b"/", 0);
        ambiguous_at::<Two>(b"abc/", 0);
    }

    #[test]
    fn two_absent() {
        absent::<Two>(b"");
        absent::<Two>(b"a/b/c");
        absent::<Two>(b"abc/xyz");
    }

    #[test]
    fn two_present_beats_ambiguous() {
        present_at::<Two>(b"//", 0);
        present_at::<Two>(b"//nope/", 0);
    }

    // ---- 5-byte delim: "ABCDE" ----
    dummy_parser!(Five, b"ABCDE");

    #[test]
    fn five_present() {
        present_at::<Five>(b"ABCDE", 0);
        present_at::<Five>(b"xxABCDEyy", 2);
        present_at::<Five>(b"AABCDE", 1);
    }

    #[test]
    fn five_ambiguous_tail() {
        ambiguous_at::<Five>(b"A", 0);
        ambiguous_at::<Five>(b"xxA", 0);
        //                        3210
        ambiguous_at::<Five>(b"xyzABCD", 3);
        //                          3210
        ambiguous_at::<Five>(b"AAA AAAAA", 3);
    }

    #[test]
    fn five_absent() {
        absent::<Five>(b"");
        absent::<Five>(b"xyz");
        absent::<Five>(b"xxxABXYZ");
        absent::<Five>(b"xxxABCDX");
    }

    #[test]
    fn five_present_beats_ambiguous() {
        present_at::<Five>(b"ABCDExxx", 0);
        present_at::<Five>(b"xxABCDE", 2);
    }
}
