use self::consts::*;

use std::{
    fs::create_dir_all,
    io::Write,
    path::{Path, PathBuf},
    process::Command,
};

use anyhow::{anyhow, bail};

use libnail::{align::Nats, alphabet::UTF8_TO_DIGITAL_AMINO};

use crate::{
    args::SearchArgs,
    io::{Database, Fasta, P7Hmm},
    util::{CommandExt, PathBufExt},
};

pub mod consts {

    pub const AMINO_DBTYPE: &[u8] = &[0, 0, 0, 0];
    pub const PROFILE_DBTYPE: &[u8] = &[2, 0, 0, 0];
    pub const GENERIC_DBTYPE: &[u8] = &[12, 0, 0, 0];
    #[cfg(not(target_os = "windows"))]
    pub const BLOSUM_62: &str = include_str!("../mat/blosum62.mat");
    #[cfg(not(target_os = "windows"))]
    pub const BLOSUM_80: &str = include_str!("../mat/blosum80.mat");
    #[cfg(target_os = "windows")]
    pub const BLOSUM_62: &str = include_str!("..\\mat\\blosum62.mat");
    #[cfg(target_os = "windows")]
    pub const BLOSUM_80: &str = include_str!("..\\mat\\blosum80.mat");
}

pub enum MmseqsScoreModel {
    Profile,
    Blosum62,
    Blosum80,
}

impl MmseqsScoreModel {
    pub fn write(&self, dir: impl AsRef<Path>) -> anyhow::Result<Option<PathBuf>> {
        let dir = dir.as_ref();

        if !dir.is_dir() {
            bail!("path: {} is not a directory", dir.to_string_lossy());
        }

        let (mat_str, mat_path) = match self {
            MmseqsScoreModel::Profile => return Ok(None),
            MmseqsScoreModel::Blosum62 => (BLOSUM_62, dir.join("blosum62.out")),
            MmseqsScoreModel::Blosum80 => (BLOSUM_80, dir.join("blosum80.out")),
        };

        mat_path.open(true)?.write_all(mat_str.as_bytes())?;

        Ok(Some(mat_path))
    }
}

pub struct MmseqsDbPaths {
    pub query_db: PathBuf,
    pub target_db: PathBuf,
    pub prefilter_db: PathBuf,
    pub align_db: PathBuf,
}

impl MmseqsDbPaths {
    pub fn new(dir: impl AsRef<Path>) -> Self {
        let dir = dir.as_ref().to_path_buf();
        Self {
            query_db: dir.join("queryDB"),
            target_db: dir.join("targetDB"),
            prefilter_db: dir.join("prefilterDB"),
            align_db: dir.join("alignDB"),
        }
    }

    pub fn dir(&self) -> anyhow::Result<&Path> {
        match self.query_db.parent() {
            Some(path) => Ok(path),
            None => bail!("no parent directory found in call to MmseqsDbPaths::write()"),
        }
    }
}

pub fn write_mmseqs_sequence_database(
    sequences: &Fasta,
    path: impl AsRef<Path>,
) -> anyhow::Result<()> {
    let db_path = path.as_ref().to_owned();
    let db_name = db_path.file_name().unwrap().to_str().unwrap();
    let db_index_path = db_path.with_file_name(format!("{db_name}.index"));
    let db_dbtype_path = db_path.with_file_name(format!("{db_name}.dbtype"));

    let header_path = db_path.with_file_name(format!("{db_name}_h"));
    let header_index_path = db_path.with_file_name(format!("{db_name}_h.index"));
    let header_dbtype_path = db_path.with_file_name(format!("{db_name}_h.dbtype"));

    let dir = db_path.parent().unwrap();
    create_dir_all(dir)?;

    db_dbtype_path.open(true)?.write_all(AMINO_DBTYPE)?;
    header_dbtype_path.open(true)?.write_all(GENERIC_DBTYPE)?;

    let mut db = db_path.open(true)?;
    let mut db_index = db_index_path.open(true)?;
    let mut db_header = header_path.open(true)?;
    let mut db_header_index = header_index_path.open(true)?;

    let mut db_offset = 0usize;
    let mut header_offset = 0usize;

    for (seq_count, seq) in sequences.values().enumerate() {
        db.write_all(&seq.utf8_bytes[1..])?;
        db.write_all(&[10u8, 0u8])?;

        let db_byte_length = seq.length + 2;
        writeln!(db_index, "{}\t{}\t{}", seq_count, db_offset, db_byte_length,)?;

        writeln!(db_header, "{}", seq.name)?;
        db_header.write_all(&[0u8])?;

        let header_byte_length = seq.name.len() + 2;
        writeln!(
            db_header_index,
            "{}\t{}\t{}",
            seq_count, header_offset, header_byte_length
        )?;

        db_offset += db_byte_length;
        header_offset += header_byte_length;
    }

    Ok(())
}

pub fn write_mmseqs_profile_database(
    profiles: &P7Hmm,
    path: impl AsRef<Path>,
) -> anyhow::Result<()> {
    let db_path = path.as_ref().to_owned();
    let db_name = db_path.file_name().unwrap().to_str().unwrap();
    let db_index_path = db_path.with_file_name(format!("{db_name}.index"));
    let db_dbtype_path = db_path.with_file_name(format!("{db_name}.dbtype"));

    let header_path = db_path.with_file_name(format!("{db_name}_h"));
    let header_index_path = db_path.with_file_name(format!("{db_name}_h.index"));
    let header_dbtype_path = db_path.with_file_name(format!("{db_name}_h.dbtype"));

    let dir = db_path.parent().unwrap();
    create_dir_all(dir)?;

    db_dbtype_path.open(true)?.write_all(PROFILE_DBTYPE)?;
    header_dbtype_path.open(true)?.write_all(GENERIC_DBTYPE)?;

    let mut db = db_path.open(true)?;
    let mut db_index = db_index_path.open(true)?;
    let mut db_header = header_path.open(true)?;
    let mut db_header_index = header_index_path.open(true)?;

    let mut db_offset = 0usize;
    let mut header_offset = 0usize;

    for (prf_cnt, prf) in profiles.values().enumerate() {
        for prf_idx in 1..=prf.length {
            for byte in (0..20)
                .map(|residue| Nats(prf.match_score(residue, prf_idx)))
                .map(|nats| nats.to_bits())
                .map(|bits| bits.value())
                // multiply the bits by 8.0 just because
                // that's what they do in mmseqs
                .map(|s| (s * 8.0) as i8)
                .map(|s| s as u8)
            {
                db.write_all(&[byte])?;
            }

            let consensus_byte_digital = *UTF8_TO_DIGITAL_AMINO
                .get(&prf.consensus_seq_bytes_utf8[prf_idx])
                .unwrap();

            // query sequence byte?
            // still not sure what this is, so we'll just put the consensus
            db.write_all(&[consensus_byte_digital])?;

            // consensus sequence byte
            db.write_all(&[consensus_byte_digital])?;

            // neff value
            db.write_all(&[0u8])?;

            // something
            db.write_all(&[0u8])?;
            // something
            db.write_all(&[0u8])?;
        }

        let db_byte_length = prf.length * 25;

        writeln!(db_index, "{}\t{}\t{}", prf_cnt, db_offset, db_byte_length,)?;

        // for some reason, the header has newlines and 0-byte separators?
        writeln!(db_header, "{}", prf.name)?;
        db_header.write_all(&[0u8])?;

        // +1 for the 0 byte, +1 for the newline
        let header_byte_length = prf.name.len() + 2;

        writeln!(
            db_header_index,
            "{}\t{}\t{}",
            prf_cnt, header_offset, header_byte_length
        )?;

        db_offset += db_byte_length;
        header_offset += header_byte_length;
    }

    Ok(())
}

pub fn run_mmseqs_search<P: AsRef<Path>>(
    paths: &MmseqsDbPaths,
    align_tsv: P,
    args: &SearchArgs,
    score_model: MmseqsScoreModel,
) -> anyhow::Result<()> {
    let effective_e_value = args.pipeline_args.seed_pvalue_threshold
        * args
            .expert_args
            .target_database_size
            .ok_or(anyhow!("no target database size"))? as f64;

    let _ = paths.align_db.remove();
    let _ = paths.align_db.with_extension("dbtype").remove();
    let _ = paths.align_db.with_extension("index").remove();

    let score_mat_path = score_model.write(paths.dir()?)?;

    {
        let mut prefilter = Command::new("mmseqs");

        prefilter
            .arg("prefilter")
            .arg(&paths.query_db)
            .arg(&paths.target_db)
            .arg(&paths.prefilter_db)
            .args(["--threads", &args.num_threads.to_string()])
            .args(["-k", &args.mmseqs_args.k.to_string()])
            .args(["-s", &args.mmseqs_args.s.to_string()])
            .args(["--max-seqs", &args.mmseqs_args.max_seqs.to_string()]);

        if let Some(v) = args.mmseqs_args.comp_bias_corr {
            prefilter.args(["--comp-bias-corr", &v.to_string()]);
        }

        if let Some(ref path) = score_mat_path {
            prefilter.arg("--sub-mat");
            prefilter.arg(path);
        };

        prefilter.run()?;
    }

    {
        let mut align = Command::new("mmseqs");
        align
            .arg("align")
            .arg(&paths.query_db)
            .arg(&paths.target_db)
            .arg(&paths.prefilter_db)
            .arg(&paths.align_db)
            .args(["--threads", &args.num_threads.to_string()])
            .args(["-e", &effective_e_value.to_string()])
            // the '-a' argument enables alignment backtraces in mmseqs2
            // it is required to get start positions for alignments
            .args(["-a", "1"]);

        if let Some(v) = args.mmseqs_args.comp_bias_corr {
            align.args(["--comp-bias-corr", &v.to_string()]);
        }

        if let Some(path) = score_mat_path {
            align.arg("--sub-mat");
            align.arg(path);
        };

        align.run()?;
    }

    Command::new("mmseqs")
        .arg("convertalis")
        .arg(&paths.query_db)
        .arg(&paths.target_db)
        .arg(&paths.align_db)
        .arg(align_tsv.as_ref())
        .args(["--threads", &args.num_threads.to_string()])
        .args([
            "--format-output",
            "qheader,theader,qstart,qend,tstart,tend,bits",
        ])
        .run()?;

    Ok(())
}
