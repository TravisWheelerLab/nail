use std::io::{Read, Seek, SeekFrom};

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
