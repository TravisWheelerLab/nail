use std::fmt::{Display, Formatter};
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

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
