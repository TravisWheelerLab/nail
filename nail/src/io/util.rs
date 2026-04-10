use std::io::{Read, Seek, SeekFrom};

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
