use anyhow::{anyhow, bail, Context, Result};
use regex::Regex;
use std::{
    fmt::{Display, Formatter},
    fs::File,
    io::{BufRead, BufReader, Read},
    path::Path,
};
use thiserror::Error;

#[derive(Default, Debug, Clone, PartialEq)]
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

#[derive(Error, Debug, PartialEq)]
#[error("can't guess file format of: {path}")]
pub struct UnrecognizedFileFormatError {
    path: String,
}

pub fn guess_query_format_from_query_file(
    query_path: &impl AsRef<Path>,
) -> Result<FileFormat> {
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

pub fn my_guess_query_format_from_query_file(
    query_path: &impl AsRef<Path>,
) -> Result<FileFormat> {
    let file = File::open(query_path).map_err(|e| {
        anyhow!("{}: {e}", query_path.as_ref().to_string_lossy())
    })?;

    let mut reader = BufReader::new(file);
    let mut first_line = String::new();
    reader.read_line(&mut first_line)?;

    let fasta_re = Regex::new("^>").unwrap();
    let sto_re = Regex::new("^# STOCKHOLM").unwrap();
    let hmmer_re = Regex::new("^HMMER").unwrap();

    if fasta_re.is_match(&first_line) {
        Ok(FileFormat::Fasta)
    } else if sto_re.is_match(&first_line) {
        Ok(FileFormat::Stockholm)
    } else if hmmer_re.is_match(&first_line) {
        Ok(FileFormat::Hmm)
    } else {
        bail!(format!(
            r#"Unknown format: "{}""#,
            query_path.as_ref().to_string_lossy().to_string()
        ))
    }
}

#[derive(Error, Debug)]
#[error("unsupported MMseqs2 database type: {code}")]
pub struct UnsupportedMmseqsDbError {
    code: u8,
}

// Not used?
pub fn _read_query_format_from_mmseqs_query_db(
    query_db_path: &impl AsRef<Path>,
) -> Result<FileFormat> {
    let mut file = File::open(query_db_path)?;
    let mut dbtype_buf: Vec<u8> = vec![];
    file.read_to_end(&mut dbtype_buf)?;
    //  from mmseqs2: commons/parameters.h
    //      DBTYPE_AMINO_ACIDS = 0;
    //      DBTYPE_NUCLEOTIDES = 1;
    //      DBTYPE_HMM_PROFILE = 2;
    //      //DBTYPE_PROFILE_STATE_SEQ = 3;
    //      //DBTYPE_PROFILE_STATE_PROFILE = 4;
    //      DBTYPE_ALIGNMENT_RES = 5;
    //      DBTYPE_CLUSTER_RES = 6;
    //      DBTYPE_PREFILTER_RES = 7;
    //      DBTYPE_TAXONOMICAL_RESULT = 8;
    //      DBTYPE_INDEX_DB = 9;
    //      DBTYPE_CA3M_DB = 10;
    //      DBTYPE_MSA_DB = 11;
    //      DBTYPE_GENERIC_DB = 12;
    //      DBTYPE_OMIT_FILE = 13;
    //      DBTYPE_PREFILTER_REV_RES = 14;
    //      DBTYPE_OFFSETDB = 15;
    //      DBTYPE_DIRECTORY = 16; // needed for verification
    //      DBTYPE_FLATFILE = 17; // needed for verification
    //      DBTYPE_SEQTAXDB = 18; // needed for verification
    //      DBTYPE_STDIN = 19; // needed for verification
    //      DBTYPE_URI = 20; // needed for verification
    match dbtype_buf[0] {
        0u8 => Ok(FileFormat::Fasta),
        2u8 => Ok(FileFormat::Stockholm),
        _ => Err(UnsupportedMmseqsDbError {
            code: dbtype_buf[0],
        }
        .into()),
    }
}

#[cfg(test)]
mod tests {
    use super::{my_guess_query_format_from_query_file, FileFormat};
    use pretty_assertions::assert_eq;
    use std::path::PathBuf;

    #[test]
    fn test_guess_query_format_from_query_file() {
        let res =
            my_guess_query_format_from_query_file(&PathBuf::from("./bad"));
        assert!(res.is_err());
        assert_eq!(
            res.unwrap_err().to_string(),
            "./bad: No such file or directory (os error 2)"
        );

        let empty = "tests/inputs/empty";
        let res =
            my_guess_query_format_from_query_file(&PathBuf::from(empty));
        assert!(res.is_err());
        let res = res.unwrap_err();
        assert_eq!(
            res.to_string(),
            r#"Unknown format: "tests/inputs/empty""#
        );

        let res = my_guess_query_format_from_query_file(&PathBuf::from(
            "tests/inputs/query.fa",
        ));
        assert!(res.is_ok());
        assert_eq!(res.unwrap(), FileFormat::Fasta);

        let res = my_guess_query_format_from_query_file(&PathBuf::from(
            "tests/inputs/query.sto",
        ));
        assert!(res.is_ok());
        assert_eq!(res.unwrap(), FileFormat::Stockholm);

        let res = my_guess_query_format_from_query_file(&PathBuf::from(
            "tests/inputs/query.hmm",
        ));
        assert!(res.is_ok());
        assert_eq!(res.unwrap(), FileFormat::Hmm);
    }
}
