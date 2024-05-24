use std::fs::File;
use std::io::BufWriter;
use std::path::PathBuf;
use std::process::Command;

use anyhow::{Context, Result};
use thiserror::Error;

#[derive(Error, Debug)]
#[error("command exited without success")]
struct CommandExitStatusError;

/// An extension trait that is intended to add a run method to the std::process::Command struct.
pub trait CommandExt {
    fn run(&mut self) -> Result<()>;
}

impl CommandExt for Command {
    fn run(&mut self) -> Result<()> {
        let output = self.output().context("failed to run command")?;

        match output.status.success() {
            true => Ok(()),
            false => {
                let stdout = std::str::from_utf8(&output.stdout)
                    .context("failed to convert sdtout to UTF8")?;
                let stderr = std::str::from_utf8(&output.stderr)
                    .context("failed to convert sdterr to UTF8")?;

                println!("command:\n{self:?}\n");
                println!("stdout:\n{stdout}\n");
                println!("stderr:\n{stderr}\n");
                Err(CommandExitStatusError.into())
            }
        }
    }
}

pub trait PathBufExt {
    fn open(&self, allow_overwrite: bool) -> anyhow::Result<BufWriter<File>>;
    fn remove(&self) -> anyhow::Result<()>;
}

impl PathBufExt for PathBuf {
    fn open(&self, allow_overwrite: bool) -> anyhow::Result<BufWriter<File>> {
        let mut file_options = File::options();

        if allow_overwrite {
            file_options.write(true).truncate(true).create(true);
        } else {
            file_options.write(true).create_new(true);
        };

        let file = file_options
            .open(self)
            .context(format!("failed to create file: {}", self.to_string_lossy()))?;

        Ok(BufWriter::new(file))
    }

    fn remove(&self) -> anyhow::Result<()> {
        std::fs::remove_file(self)?;
        Ok(())
    }
}
