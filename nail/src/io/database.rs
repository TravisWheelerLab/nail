use std::{
    fs::File,
    io::{BufWriter, Read, Seek, SeekFrom, Write},
    path::{Path, PathBuf},
    sync::{Arc, Mutex},
};

use crate::io::{ReadSeekExt, ReadState};

use anyhow::{anyhow, bail, Context};
use indexmap::IndexMap;
use rayon::iter::{IntoParallelRefMutIterator, ParallelIterator};

use std::fmt::Display;

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

pub trait RecordParser: Clone + Send + Sync + 'static {
    const DELIM: &'static [u8];
    const DELIM_HEAD: u8 = Self::DELIM[0];
    const DELIM_TAIL: &'static [u8] = tail(Self::DELIM);
    const DELIM_LEN: usize = Self::DELIM.len();
    const DELIM_TAIL_LEN: usize = Self::DELIM_TAIL.len();
    const DELIM_TYPE: Delimiter;
    type Record: Display + Send + Sync;

    fn new(start_pos: u64) -> Self;
    fn offset(&mut self, line: &[u8], start: u64) -> Option<(String, Offset)>;
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

#[derive(Default, Clone, Debug, PartialEq)]
pub struct Offset {
    pub start: u64,
    pub n_bytes: usize,
}

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
            self.key, self.first.start, self.first.n_bytes, self.second.start, self.second.n_bytes
        )
    }
}

impl std::error::Error for IndexCollisionError {}

// rayon expects the the iterators it
// uses to implement DoubleEndedIterator
impl<'a, R, K> DoubleEndedIterator for DatabaseIter<'a, R, K>
where
    R: RecordParser,
    K: DoubleEndedIterator<Item = &'a str>,
{
    fn next_back(&mut self) -> Option<Self::Item> {
        self.keys.next_back().and_then(|n| self.inner.get(n))
    }
}

impl<'a, R, K> ExactSizeIterator for DatabaseIter<'a, R, K>
where
    R: RecordParser,
    K: DoubleEndedIterator<Item = &'a str>,
{
}

pub struct DatabaseIter<'a, R, K>
where
    R: RecordParser,
    K: DoubleEndedIterator<Item = &'a str>,
{
    pub(crate) inner: Database<R>,
    pub(crate) keys: K,
}

impl<'a, R, K> Iterator for DatabaseIter<'a, R, K>
where
    R: RecordParser,
    K: DoubleEndedIterator<Item = &'a str>,
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

pub type DefaultIndex = Index<IndexMap<String, Offset>>;

pub struct Database<R>
where
    R: RecordParser,
{
    pub(crate) path: PathBuf,
    pub(crate) file: File,
    pub(crate) index: Arc<DefaultIndex>,
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
            file: File::open(&self.path).expect("failed to re-open database file"),
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
        let index = Index::build::<_, R, 15, 1>(path.as_ref())?;
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

    // note: this is called a Return-position impl trait, acting as an abstract return type
    //       see: https://doc.rust-lang.org/reference/types/impl-trait.html
    pub fn iter(&'_ self) -> DatabaseIter<'_, R, impl DoubleEndedIterator<Item = &str>> {
        DatabaseIter {
            inner: self.clone(),
            keys: self.index.keys(),
        }
    }

    pub fn write<W: Write>(&self, out: W) -> anyhow::Result<()> {
        let mut out = BufWriter::new(out);
        self.iter().try_for_each(|rec| {
            let rec = rec.context("itertor failed to produce a record")?;
            writeln!(out, "{rec}").context("failed to write record in call to Database::write()")
        })?;

        Ok(())
    }

    pub fn get(&mut self, key: &str) -> Option<anyhow::Result<R::Record>> {
        let offset = self.index.get(key)?;

        self.buffer.resize(offset.n_bytes, 0u8);

        if let Err(e) = self.file.seek(SeekFrom::Start(offset.start)) {
            return Some(Err(e).context("failed to seek to offset in database file"));
        }

        if let Err(e) = self.file.read_exact(&mut self.buffer) {
            return Some(Err(e).context("failed to read from database file"));
        };

        match R::parse(&self.buffer) {
            Ok(rec) => Some(Ok(rec)),
            Err(e) => Some(Err(e).context("failed to parse record")),
        }
    }
}

pub trait IndexInner {
    fn new() -> Self;
    fn len(&self) -> usize;
    // note: this is called a Return-position impl trait in trait (RPITIT)
    //       it's kind of automatically translated into some (hidden) concrete
    //       type depending on the implementer
    //
    //  see: https://rustc-dev-guide.rust-lang.org/return-position-impl-trait-in-trait.html
    fn keys(&self) -> impl DoubleEndedIterator<Item = &str>;
    fn get(&self, key: &str) -> Option<&Offset>;
    fn extend<I: IntoIterator<Item = (String, Offset)>>(
        &mut self,
        iterable: I,
    ) -> anyhow::Result<()>;
}

impl IndexInner for IndexMap<String, Offset> {
    fn new() -> Self {
        Self::new()
    }

    fn len(&self) -> usize {
        self.len()
    }

    fn keys(&self) -> impl DoubleEndedIterator<Item = &str> {
        self.keys().map(|k| k.as_str())
    }

    fn get(&self, key: &str) -> Option<&Offset> {
        self.get(key)
    }

    fn extend<I: IntoIterator<Item = (String, Offset)>>(
        &mut self,
        iterable: I,
    ) -> anyhow::Result<()> {
        let iter = iterable.into_iter();

        let (lower_len, _) = iter.size_hint();
        let reserve = if self.is_empty() {
            lower_len
        } else {
            lower_len.div_ceil(2)
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

pub struct Index<I>
where
    I: IndexInner,
{
    inner: I,
}

impl<I> Index<I>
where
    I: IndexInner + Send + Sync,
{
    pub fn len(&self) -> usize {
        self.inner.len()
    }

    pub fn keys(&self) -> impl DoubleEndedIterator<Item = &str> {
        self.inner.keys()
    }

    pub fn get(&self, key: &str) -> Option<&Offset> {
        self.inner.get(key)
    }

    pub(crate) fn build<
        P: AsRef<Path>,
        R: RecordParser,
        const BUF_EXP: usize,
        const CHUNKS: usize,
    >(
        path: P,
    ) -> anyhow::Result<Self> {
        let mut file = File::open(path.as_ref())?;

        let sz = file.metadata().unwrap().len();
        let chunk_sz = sz / CHUNKS as u64;
        let mut delim_positions: Vec<u64> = vec![];

        // 2 << (B - 1) == 2^B
        let mut buf = vec![0; 2 << (BUF_EXP - 1)];

        use ReadState::*;
        'outer: for mut start in (1..CHUNKS).map(|i| i as u64 * chunk_sz) {
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
                block_starts.extend(block_ends.iter().take(block_ends.len() - 1).map(|b| b + 1));

                (block_starts, block_ends)
            }
        };

        let mut files_and_starts: Vec<_> = block_starts
            .into_iter()
            .zip(block_ends)
            .map(|(start, end)| {
                let mut file = File::open(path.as_ref())?;
                let limit = end - start + 1;
                file.seek(SeekFrom::Start(start))?;
                Ok((file.take(limit), start))
            })
            .collect::<anyhow::Result<_>>()?;

        let index = Arc::new(Mutex::new(I::new()));

        if let Err(e) = files_and_starts
            .par_iter_mut()
            .try_for_each(|(file, start)| {
                Self::build_chunk::<_, R, 5_000, 20>(file, index.clone(), Some(*start))
            })
        {
            match e.downcast_ref::<IndexCollisionError>() {
                Some(e) => {
                    use ReadState::*;
                    let mut buf = vec![0; 2 << (BUF_EXP - 1)];

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

        Ok(Self {
            inner: Arc::into_inner(index)
                .ok_or(anyhow!("no index"))?
                .into_inner()
                .ok()
                .context("mutex poisoned")?,
        })
    }

    fn build_chunk<R, P, const REC_LIM: usize, const BUF_EXP: usize>(
        mut data_chunk: R,
        index: Arc<Mutex<I>>,
        start: Option<u64>,
    ) -> anyhow::Result<()>
    where
        R: Read + Seek,
        P: RecordParser,
    {
        let process_fn = |recs: &mut Vec<(String, Offset)>, done: bool| {
            let full = recs.len() == REC_LIM;

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
                Some(mut guard) => guard.extend(recs.drain(..)),
                None => Ok(()),
            }
        };

        use ReadState::*;
        let mut pos = start.unwrap_or(0);
        let mut parser = P::new(pos);
        let mut on_first_record = true;
        let mut buf = vec![0; 2 << (BUF_EXP - 1)];
        let mut names_and_offsets: Vec<(String, Offset)> = Vec::with_capacity(REC_LIM);

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

    macro_rules! dummy_parser {
        ($name:ident, $delim:expr) => {
            #[derive(Debug, Clone)]
            struct $name;

            impl RecordParser for $name {
                const DELIM: &'static [u8] = $delim;
                const DELIM_TYPE: Delimiter = Delimiter::Initiating;
                type Record = String;

                fn new(_: u64) -> Self {
                    $name
                }
                fn offset(&mut self, _: &[u8], _: u64) -> Option<(String, Offset)> {
                    None
                }
                fn parse(_: &[u8]) -> anyhow::Result<Self::Record> {
                    Ok("".to_string())
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
