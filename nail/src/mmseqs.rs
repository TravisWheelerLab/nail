use self::consts::*;

use std::{
    collections::HashMap,
    fs::{create_dir_all, File},
    io::{BufRead, BufReader, Write},
    path::{Path, PathBuf},
    process::Command,
};

use anyhow::Context;
use clap::Args;

use libnail::{
    align::{structs::Seed, Nats},
    alphabet::UTF8_TO_DIGITAL_AMINO,
    structs::{Profile, Sequence},
};

use crate::{
    pipeline::SeedMap,
    util::{CommandExt, PathBufExt},
};

pub mod consts {
    pub const AMINO_DBTYPE: &[u8] = &[0, 0, 0, 0];
    pub const PROFILE_DBTYPE: &[u8] = &[2, 0, 0, 0];
    pub const GENERIC_DBTYPE: &[u8] = &[12, 0, 0, 0];
}

#[derive(Args, Debug, Clone, Default)]
pub struct MmseqsArgs {
    /// The directory where intermediate files will be placed
    #[arg(long = "prep", value_name = "<PATH>", default_value = "prep/")]
    pub prep_dir: PathBuf,

    /// MMseqs2 prefilter: k-mer length (0: automatically set to optimum)
    #[arg(long = "mmseqs-k", default_value_t = 0usize)]
    pub k: usize,

    /// MMseqs2 prefilter: k-mer threshold for generating similar k-mer lists
    #[arg(long = "mmseqs-k-score", default_value_t = 80usize)]
    pub k_score: usize,

    /// MMseqs2 prefilter: Accept only matches with ungapped alignment score above threshold
    #[arg(long = "mmseqs-min-ungapped_score", default_value_t = 15usize)]
    pub min_ungapped_score: usize,

    /// MMseqs2 prefilter: Maximum results per query sequence allowed to pass the prefilter
    #[arg(long = "mmseqs-max-seqs", default_value_t = 1000usize)]
    pub max_seqs: usize,

    /// MMseqs2 align: Include matches below this P-value as seeds.
    ///
    /// Note: the MMseqs2 align tool only allows thresholding by E-value, so the P-value supplied
    /// here is multiplied by the size of the target database (i.e. number of sequences) to achieve
    /// an E-value threshold that is effectively the same as the chosen P-value threshold.
    #[arg(long = "mmseqs-pvalue-threshold", default_value_t = 0.01f64)]
    pub pvalue_threshold: f64,
}

pub fn write_mmseqs_sequence_database(
    sequences: &[Sequence],
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

    for (seq_count, seq) in sequences.iter().enumerate() {
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
    profiles: &[impl AsRef<Profile>],
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

    for (profile_count, profile) in profiles.iter().map(|p| p.as_ref()).enumerate() {
        for profile_idx in 1..=profile.length {
            for byte in (0..20)
                .map(|residue| Nats(profile.match_score(residue, profile_idx)))
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
                .get(&profile.consensus_sequence_bytes_utf8[profile_idx])
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

        let db_byte_length = profile.length * 25;

        writeln!(
            db_index,
            "{}\t{}\t{}",
            profile_count, db_offset, db_byte_length,
        )?;

        // for some reason, the header has newlines and 0-byte separators?
        writeln!(db_header, "{}", profile.name)?;
        db_header.write_all(&[0u8])?;

        // +1 for the 0 byte, +1 for the newline
        let header_byte_length = profile.name.len() + 2;

        writeln!(
            db_header_index,
            "{}\t{}\t{}",
            profile_count, header_offset, header_byte_length
        )?;

        db_offset += db_byte_length;
        header_offset += header_byte_length;
    }

    Ok(())
}

pub fn run_mmseqs_search(
    query_db: impl AsRef<Path>,
    target_db: impl AsRef<Path>,
    prefilter_db: impl AsRef<Path>,
    align_db: impl AsRef<Path>,
    align_tsv: impl AsRef<Path>,
    num_targets: usize,
    num_threads: usize,
    args: &MmseqsArgs,
) -> anyhow::Result<()> {
    let effective_e_value = args.pvalue_threshold * num_targets as f64;

    let query_db = query_db.as_ref();
    let target_db = target_db.as_ref();
    let prefilter_db = prefilter_db.as_ref();
    let align_db = align_db.as_ref().to_path_buf();
    let align_tsv = align_tsv.as_ref();

    let _ = align_db.remove();
    let _ = align_db.with_extension("dbtype").remove();
    let _ = align_db.with_extension("index").remove();

    Command::new("mmseqs")
        .arg("prefilter")
        .arg(query_db)
        .arg(target_db)
        .arg(prefilter_db)
        .args(["--threads", &num_threads.to_string()])
        .args(["-k", &args.k.to_string()])
        .args(["--k-score", &args.k_score.to_string()])
        .args(["--min-ungapped-score", &args.min_ungapped_score.to_string()])
        .args(["--max-seqs", &args.max_seqs.to_string()])
        .run()?;

    Command::new("mmseqs")
        .arg("align")
        .arg(query_db)
        .arg(target_db)
        .arg(prefilter_db)
        .arg(&align_db)
        .args(["--threads", &num_threads.to_string()])
        .args(["-e", &effective_e_value.to_string()])
        // the '-a' argument enables alignment backtraces in mmseqs2
        // it is required to get start positions for alignments
        .args(["-a", "1"])
        .run()?;

    Command::new("mmseqs")
        .arg("convertalis")
        .arg(query_db)
        .arg(target_db)
        .arg(align_db)
        .arg(align_tsv)
        .args(["--threads", &num_threads.to_string()])
        .args([
            "--format-output",
            "qheader,theader,qstart,qend,tstart,tend,evalue",
        ])
        .run()?;

    Ok(())
}

pub fn seeds_from_mmseqs_align_tsv(path: impl AsRef<Path>) -> anyhow::Result<SeedMap> {
    let path = path.as_ref();

    let mut seed_map: SeedMap = HashMap::new();

    let mmseqs_align_file = File::open(path).context(format!(
        "couldn't open mmseqs align file at: {}",
        path.to_string_lossy()
    ))?;

    let align_reader = BufReader::new(mmseqs_align_file);

    for line in align_reader.lines().map_while(Result::ok) {
        let line_tokens: Vec<&str> = line.split('\t').collect();

        let target_header = line_tokens[1];
        let target_header_tokens: Vec<&str> = target_header.split_whitespace().collect();
        let target_name = target_header_tokens[0].to_string();
        let target_start = line_tokens[4].parse::<usize>()?;
        let target_end = line_tokens[5].parse::<usize>()?;

        let query_header = line_tokens[0];
        let query_header_tokens: Vec<&str> = query_header.split_whitespace().collect();
        let profile_name = query_header_tokens[0].to_string();
        let profile_start = line_tokens[2].parse::<usize>()?;
        let profile_end = line_tokens[3].parse::<usize>()?;

        let profile_map = seed_map.entry(profile_name).or_default();
        profile_map.insert(
            target_name,
            Seed {
                target_start,
                target_end,
                profile_start,
                profile_end,
            },
        );
    }

    Ok(seed_map)
}
