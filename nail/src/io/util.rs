use std::io::{self, Read, Seek, SeekFrom};

use anyhow::anyhow;

#[allow(dead_code)]
pub(crate) trait ByteBufferExt {
    /// Get the index of the first non-whitespace byte
    fn first_non_whitespace_pos(&self) -> Option<usize>;
    /// Grab the string up to the next whitespace starting from pos
    fn word_from(&self, pos: usize) -> anyhow::Result<&str>;
    /// Grab the first word, skipping leading whitespace
    fn first_word(&self) -> anyhow::Result<&str>;
    /// Grab the string between start and end
    fn str(&self, start: usize, end: usize) -> anyhow::Result<&str>;
    /// Interpret the buffer as &str
    fn as_str(&self) -> anyhow::Result<&str>;
}

impl ByteBufferExt for &[u8] {
    fn word_from(&self, pos: usize) -> anyhow::Result<&str> {
        let end = self[pos..]
            .iter()
            .position(|b| b.is_ascii_whitespace())
            .map(|n| pos + n)
            .unwrap_or(self.len());

        Ok(std::str::from_utf8(&self[pos..end])?)
    }

    fn first_non_whitespace_pos(&self) -> Option<usize> {
        self.iter().position(|b| !b.is_ascii_whitespace())
    }

    fn first_word(&self) -> anyhow::Result<&str> {
        self.word_from(
            self.first_non_whitespace_pos()
                .ok_or(anyhow!("no words found in buffer"))?,
        )
    }

    fn str(&self, start: usize, end: usize) -> anyhow::Result<&str> {
        Ok(std::str::from_utf8(&self[start..=end])?)
    }

    fn as_str(&self) -> anyhow::Result<&str> {
        Ok(std::str::from_utf8(&self[0..self.len()])?)
    }
}

#[derive(Debug)]
pub enum ReadState {
    Reading(usize),
    Final(usize),
    Done,
}

pub trait ReadSeekExt: Read + Seek {
    /// Calls read on `inner` and returns a `ReadState`
    fn read_with_state(&mut self, buf: &mut [u8]) -> anyhow::Result<ReadState> {
        let n_bytes_read = self.read(buf)?;
        let pos = self.stream_position()?;
        let end = self.seek(SeekFrom::End(0))?;
        self.seek(SeekFrom::Start(pos))?;

        let finished = pos >= end;

        Ok(match (n_bytes_read, finished) {
            // read no bytes, at the end
            (0, true) => ReadState::Done,
            // read no bytes, not at the end
            (0, false) => ReadState::Reading(0),
            // read some bytes, at the end
            (n, true) => ReadState::Final(n),
            // read some bytes, not at the end
            (n, false) => ReadState::Reading(n),
        })
    }
}

impl<T: Read + Seek> ReadSeekExt for T {}

pub struct SeekableTake<T> {
    inner: T,
    /// The offset into `inner` at which the `SeekableTake` starts
    start: u64,
    /// The maximum number of bytes we are 'taking' after `start`
    limit: u64,
    /// The current relative position in the buffer (starts at 0)
    pos: u64,
}

impl<T: Read + Seek> SeekableTake<T> {
    pub fn new(mut inner: T, start: u64, limit: u64) -> io::Result<Self> {
        let end = inner.seek(SeekFrom::End(0))?;
        let start = start.min(end);
        let limit = limit.min(end.saturating_sub(start));
        inner.seek(SeekFrom::Start(start))?;

        Ok(Self {
            inner,
            start,
            limit,
            pos: 0,
        })
    }

    fn remaining(&self) -> u64 {
        self.limit.saturating_sub(self.pos)
    }
}

impl<T: Read + Seek> Read for SeekableTake<T> {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        // read at most self.remaining() bytes
        let n_bytes_to_read = (buf.len() as u64).min(self.remaining()) as usize;
        let n_bytes_read = self.inner.read(&mut buf[..n_bytes_to_read])?;
        self.pos += n_bytes_read as u64;
        Ok(n_bytes_read)
    }
}

impl<T: Read + Seek> Seek for SeekableTake<T> {
    fn seek(&mut self, pos: SeekFrom) -> io::Result<u64> {
        // map the pos relative to start
        let new_pos = match pos {
            SeekFrom::Start(n) => Some(n),
            SeekFrom::End(n) => (self.limit as i64).checked_add(n).map(|n| n as u64),
            SeekFrom::Current(n) => (self.pos as i64).checked_add(n).map(|n| n as u64),
        }
        .ok_or_else(|| {
            io::Error::new(io::ErrorKind::InvalidInput, "invalid seek in SeekableTake")
        })?;

        if new_pos > self.limit {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "seek past limit",
            ));
        }

        self.inner.seek(SeekFrom::Start(self.start + new_pos))?;
        self.pos = new_pos;
        Ok(new_pos)
    }
}
