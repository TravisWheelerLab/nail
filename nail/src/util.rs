use std::fmt::{Display, Formatter};
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Read, Seek, SeekFrom};
use std::path::Path;
use std::process::Command;
use std::time::Instant;

use anyhow::{bail, Context};
use thiserror::Error;

#[allow(dead_code)]
pub mod term {
    pub const RESET: &str = "\x1b[0m";

    pub const RED: &str = "\x1b[31m";
    pub const GREEN: &str = "\x1b[32m";
    pub const YELLOW: &str = "\x1b[33m";
    pub const BLUE: &str = "\x1b[34m";
    pub const MAGENTA: &str = "\x1b[35m";
    pub const CYAN: &str = "\x1b[36m";

    pub const BRIGHT_RED: &str = "\x1b[91m";
    pub const BRIGHT_GREEN: &str = "\x1b[92m";
    pub const BRIGHT_YELLOW: &str = "\x1b[93m";
    pub const BRIGHT_BLUE: &str = "\x1b[94m";
    pub const BRIGHT_MAGENTA: &str = "\x1b[95m";
    pub const BRIGHT_CYAN: &str = "\x1b[96m";
}

#[allow(dead_code)]
pub fn burn(dur: std::time::Duration) {
    let start = Instant::now();
    while Instant::now() - start < dur {
        std::hint::spin_loop();
    }
}

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

#[allow(dead_code)]
pub trait FileExt {
    fn byte(&mut self, n: u64) -> u8;
}

impl FileExt for File {
    fn byte(&mut self, n: u64) -> u8 {
        self.seek(SeekFrom::Start(n)).expect("failed to seek");
        let mut byte = [0u8; 1];
        self.read_exact(&mut byte).expect("failed to read");
        byte[0]
    }
}

pub trait PathExt: AsRef<Path> {
    fn check_open(&self, allow_overwrite: bool) -> anyhow::Result<()>;
    fn open(&self, allow_overwrite: bool) -> anyhow::Result<BufWriter<File>>;
    fn create_dir(&self) -> anyhow::Result<()>;
}

impl<P: AsRef<Path>> PathExt for P {
    fn check_open(&self, allow_overwrite: bool) -> anyhow::Result<()> {
        let path = self.as_ref();

        let exists = std::fs::exists(path)
            .with_context(|| format!("failed to check existence of {path:?}"))?;

        match (exists, allow_overwrite) {
            (true, true) => {
                let f = File::options().write(true).open(path)?;
                let len = f.metadata()?.len();
                // checks if we can truncate without actually truncating
                f.set_len(len)?;
                Ok(())
            }
            (true, false) => bail!("file: {path:?} already exists"),
            (false, true) | (false, false) => {
                let f = File::options()
                    .write(true)
                    .create_new(true)
                    .open(path)
                    .with_context(|| format!("failed to open-check {path:?}"))?;
                drop(f);
                std::fs::remove_file(path)
                    .with_context(|| format!("failed to remove-check {path:?}"))?;
                Ok(())
            }
        }
    }

    fn open(&self, allow_overwrite: bool) -> anyhow::Result<BufWriter<File>> {
        let path = self.as_ref();
        let mut opts = File::options();

        if allow_overwrite {
            opts.write(true).truncate(true).create(true);
        } else {
            opts.write(true).create_new(true);
        };

        let file = opts
            .open(path)
            .with_context(|| format!("failed to create file: {path:?}"))?;

        Ok(BufWriter::new(file))
    }

    fn create_dir(&self) -> anyhow::Result<()> {
        let path = self.as_ref();
        if !path
            .try_exists()
            .with_context(|| format!("failed to check existence of: {path:?}"))?
        {
            std::fs::create_dir(path)
                .with_context(|| format!("failed to create directory: {path:?}"))?;
        }

        Ok(())
    }
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
