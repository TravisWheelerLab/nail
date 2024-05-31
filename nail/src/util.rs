use std::fmt::{Display, Formatter};
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter};
use std::path::{Path, PathBuf};
use std::process::Command;

use anyhow::Context;
use thiserror::Error;

#[derive(Default, Debug, Clone)]
pub enum FileFormat {
    Fasta,
    Stockholm,
    Hmm,
    #[default]
    Unset,
}

impl Display for FileFormat {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        match self {
            FileFormat::Fasta => write!(f, "Fasta"),
            FileFormat::Stockholm => write!(f, "Stockholm"),
            FileFormat::Hmm => write!(f, "HMM"),
            FileFormat::Unset => write!(f, "Unset"),
        }
    }
}

#[derive(Error, Debug)]
#[error("can't guess file format of: {path}")]
pub struct UnrecognizedFileFormatError {
    path: String,
}

pub fn guess_query_format_from_query_file(
    query_path: &impl AsRef<Path>,
) -> anyhow::Result<FileFormat> {
    let file = File::open(query_path).context(format!(
        "failed to open query file: {}",
        query_path.as_ref().to_string_lossy()
    ))?;

    let mut reader = BufReader::new(file);
    let mut first_line = String::new();
    reader.read_line(&mut first_line)?;

    if &first_line[0..1] == ">" {
        Ok(FileFormat::Fasta)
    } else if &first_line[0..11] == "# STOCKHOLM" {
        Ok(FileFormat::Stockholm)
    } else if &first_line[0..5] == "HMMER" {
        Ok(FileFormat::Hmm)
    } else {
        Err(UnrecognizedFileFormatError {
            path: query_path.as_ref().to_string_lossy().to_string(),
        }
        .into())
    }
}

#[derive(Error, Debug)]
#[error("command exited without success")]
struct CommandExitStatusError;

/// An extension trait that is intended to add a run method to the std::process::Command struct.
pub trait CommandExt {
    fn run(&mut self) -> anyhow::Result<()>;
}

impl CommandExt for Command {
    fn run(&mut self) -> anyhow::Result<()> {
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

pub fn check_hmmer_installed() -> anyhow::Result<()> {
    Command::new("hmmbuild")
        .arg("-h")
        .run()
        .context("hmmbuild does not appear to be in the system path")
}

pub fn check_mmseqs_installed() -> anyhow::Result<()> {
    Command::new("mmseqs")
        .arg("-h")
        .run()
        .context("mmseqs2 does not appear to be in the system path")
}

pub fn set_threads(num_threads: usize) -> anyhow::Result<()> {
    rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build_global()
        .context("failed to build rayon global threadpool")
}
