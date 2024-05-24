use crate::args::{guess_query_format_from_query_file, FileFormat};
use crate::extension_traits::{CommandExt, PathBufExt};

use libnail::align::structs::Seed;
use libnail::align::{needleman_wunsch, Nats, SimpleTraceStep};
use libnail::alphabet::UTF8_TO_DIGITAL_AMINO;
use libnail::structs::hmm::parse_hmms_from_p7hmm_file;
use libnail::structs::{Profile, Sequence};
use std::collections::HashMap;
use std::fs::{create_dir_all, File};
use std::io::Write;
use std::io::{BufRead, BufReader, Read, Seek, SeekFrom};
use std::path::PathBuf;
use std::process::Command;

use crate::cli::CommonArgs;
use anyhow::Context;
use clap::Args;
use thiserror::Error;

pub type SeedMap = HashMap<String, HashMap<String, Seed>>;

pub fn p7_to_mmseqs_profile(profiles: &[Profile], args: &PrepDirArgs) -> anyhow::Result<()> {
    create_dir_all(&args.path)?;

    args.mmseqs_query_h_dbtype_path()
        .open(true)?
        .write_all(&[12, 0, 0, 0])?;

    args.mmseqs_query_dbtype_path()
        .open(true)?
        .write_all(&[2, 0, 0, 0])?;

    let mut query_db = args.mmseqs_query_db_path().open(true)?;
    let mut query_db_index = args.mmseqs_query_db_index_path().open(true)?;
    let mut query_db_header = args.mmseqs_query_db_h_path().open(true)?;
    let mut query_db_header_index = args.mmseqs_query_db_h_index_path().open(true)?;

    let mut query_offset = 0usize;
    let mut header_offset = 0usize;

    for (profile_count, profile) in profiles.iter().enumerate() {
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
                query_db.write_all(&[byte])?;
            }

            let consensus_byte_digital = *UTF8_TO_DIGITAL_AMINO
                .get(&profile.consensus_sequence_bytes_utf8[profile_idx])
                .unwrap();

            // query sequence byte?
            // still not sure what this is, so we'll just put the consensus
            query_db.write_all(&[consensus_byte_digital])?;

            // consensus sequence byte
            query_db.write_all(&[consensus_byte_digital])?;

            // neff value
            query_db.write_all(&[0u8])?;

            // something
            query_db.write_all(&[0u8])?;
            // something
            query_db.write_all(&[0u8])?;
        }

        let query_byte_length = profile.length * 25;

        writeln!(
            query_db_index,
            "{}\t{}\t{}",
            profile_count, query_offset, query_byte_length,
        )?;

        // for some reason, the header has newlines and 0-byte separators?
        writeln!(query_db_header, "{}", profile.name)?;
        query_db_header.write_all(&[0u8])?;

        // +1 for the 0 byte, +1 for the newline
        let header_byte_length = profile.name.len() + 2;

        writeln!(
            query_db_header_index,
            "{}\t{}\t{}",
            profile_count, header_offset, header_byte_length
        )?;

        query_offset += query_byte_length;
        header_offset += header_byte_length;
    }

    Ok(())
}

#[derive(Args, Debug, Clone)]
pub struct MmseqsArgs {
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

#[derive(Args, Debug, Clone, Default)]
pub struct PrepDirArgs {
    /// Where MMseqs2 intermediate files are placed
    #[arg(
        short = 'p',
        long = "prep",
        default_value = "./prep/",
        value_name = "PATH"
    )]
    pub path: PathBuf,
}

impl PrepDirArgs {
    /// Produce a path to the query P7 HMM file
    pub fn prep_query_hmm_path(&self) -> PathBuf {
        self.path.join("query.hmm")
    }
    /// Produce a path to the MMseqs2 MSA database.
    ///
    /// This is created if a stockholm file query is provided.
    pub fn mmseqs_msa_db_path(&self) -> PathBuf {
        self.path.join("msaDB")
    }
    /// Produce a path to the MMseqs2 query database.
    ///
    /// If a fasta target was provided, this will be a sequence database.
    /// If a stockholm target was provided, this will be a profile database.
    pub fn mmseqs_query_db_path(&self) -> PathBuf {
        self.path.join("queryDB")
    }
    /// Produce a path to the MMseqs2 query database dbtype file.
    ///
    /// This file holds a byte that describes the original query file format.
    pub fn mmseqs_query_dbtype_path(&self) -> PathBuf {
        self.path.join("queryDB.dbtype")
    }
    /// Produce a path to the MMseqs2 query database index
    pub fn mmseqs_query_db_index_path(&self) -> PathBuf {
        self.path.join("queryDB.index")
    }
    /// Produce a path to the MMseqs2 query database header file
    pub fn mmseqs_query_db_h_path(&self) -> PathBuf {
        self.path.join("queryDB_h")
    }
    /// Produce a path to the MMseqs2 query database header file index
    pub fn mmseqs_query_db_h_index_path(&self) -> PathBuf {
        self.path.join("queryDB_h.index")
    }
    /// Produce a path to the MMseqs2 query database header dbtype file.
    ///
    /// This file holds a byte that describes the original query file format.
    pub fn mmseqs_query_h_dbtype_path(&self) -> PathBuf {
        self.path.join("queryDB_h.dbtype")
    }
    /// Produce a path to the MMseqs2 target database.
    ///
    /// This will always be a sequence database.
    pub fn mmseqs_target_db_path(&self) -> PathBuf {
        self.path.join("targetDB")
    }
    /// Produce a path to the MMseqs2 prefilter database.
    ///
    /// This is the result of running `mmseqs prefilter` on the query and target databases.
    pub fn mmseqs_prefilter_db_path(&self) -> PathBuf {
        self.path.join("prefilterDB")
    }
    /// Produce a path to the MMseqs2 alignment database.
    ///
    /// This is the result of running `mmseqs align` on the query, target, and prefilter databases.
    pub fn mmseqs_align_db_path(&self) -> PathBuf {
        self.path.join("alignDB")
    }
    /// Produce a path to the MMseqs2 alignment output.
    ///
    /// This is the result of running `mmseqs convertalis` on the query, target, and align databases.
    pub fn mmseqs_align_tsv_path(&self) -> PathBuf {
        self.path.join("align.tsv")
    }
}

#[derive(Args, Debug, Clone)]
pub struct SeedArgs {
    /// The location of files prepared with nail prep
    // NOTE: this arg is here so that prep_dir_path can be a positional argument.
    //       the value assigned to prep_dir_path needs to be turned into a
    //       PrepDirArgs struct before the seed function can be run
    #[arg(value_name = "PATH")]
    pub prep_dir_path: PathBuf,

    /// Where to place the seeds output file
    #[arg(short, long, default_value = "seeds.json")]
    pub seeds_path: PathBuf,

    /// Query file
    #[arg(value_name = "QUERY.[fasta:sto]")]
    pub query_path: PathBuf,

    /// Target file
    #[arg(value_name = "TARGET.fasta")]
    pub target_path: PathBuf,

    /// Provides paths to the prep dir and the files placed in it
    // TODO: error if this is not properly initialized
    #[clap(skip)]
    pub prep_dir: PrepDirArgs,

    /// Arguments that are passed to MMseqs2
    #[command(flatten)]
    pub mmseqs_args: MmseqsArgs,

    /// Arguments that are common across all nail subcommands
    #[command(flatten)]
    pub common_args: CommonArgs,
}

fn run_mmseqs_search(args: &SeedArgs) -> anyhow::Result<()> {
    create_dir_all(&args.prep_dir.path)?;

    let query_format = guess_query_format_from_query_file(&args.query_path)?;
    match query_format {
        FileFormat::Fasta => {
            Command::new("mmseqs")
                .arg("createdb")
                .arg(&args.query_path)
                .arg(&args.prep_dir.mmseqs_query_db_path())
                .run()?;

            // NOTE: we no longer prebuild HMMs from fasta queries
        }
        FileFormat::Hmm => {
            let profiles = parse_hmms_from_p7hmm_file(&args.query_path)
                .context("failed to read query hmm")?
                .iter()
                .map(Profile::new)
                .collect::<Vec<_>>();
            p7_to_mmseqs_profile(&profiles, &args.prep_dir)?;
        }
        FileFormat::Stockholm => {
            Command::new("mmseqs")
                .arg("convertmsa")
                // this flag should force the profile to be labeled
                // using the actual name rather than the accession number
                .args(["--identifier-field", "0"])
                .arg(&args.query_path)
                .arg(&args.prep_dir.mmseqs_msa_db_path())
                .run()?;

            Command::new("mmseqs")
                .arg("msa2profile")
                .arg(&args.prep_dir.mmseqs_msa_db_path())
                .arg(&args.prep_dir.mmseqs_query_db_path())
                .args(["--threads", &args.common_args.num_threads.to_string()])
                // --match-mode INT       0: Columns that have a residue in the first sequence are kept,
                //                        1: columns that have a residue in --match-ratio of all sequences
                //                           are kept [0]
                .args(["--match-mode", "1"])
                .run()?;

            Command::new("hmmbuild")
                .args(["--cpu", &args.common_args.num_threads.to_string()])
                .arg(&args.prep_dir.prep_query_hmm_path())
                .arg(&args.query_path)
                .run()?;
        }
        ref format => {
            return Err(InvalidFileFormatError {
                format: format.clone(),
            })
            .context("invalid query file format in nail prep")
        }
    }

    Command::new("mmseqs")
        .arg("createdb")
        .arg(&args.target_path)
        .arg(&args.prep_dir.mmseqs_target_db_path())
        .run()?;

    // count the number of lines in the target database so we can compute
    // the effective E-value that is equivalent to the chosen P-value
    let target_db_file = BufReader::new(File::open(args.prep_dir.mmseqs_target_db_path())?);
    let num_targets = target_db_file.lines().count() as f64;
    let effective_e_value = args.mmseqs_args.pvalue_threshold * num_targets;

    let _ = args.prep_dir.mmseqs_align_db_path().remove();

    let _ = args
        .prep_dir
        .mmseqs_align_db_path()
        .with_extension("dbtype")
        .remove();

    let _ = args
        .prep_dir
        .mmseqs_align_db_path()
        .with_extension("index")
        .remove();

    Command::new("mmseqs")
        .arg("search")
        .arg(&args.prep_dir.mmseqs_query_db_path())
        .arg(&args.prep_dir.mmseqs_target_db_path())
        .arg(&args.prep_dir.mmseqs_align_db_path())
        .arg(&args.prep_dir.path)
        .args(["--threads", &args.common_args.num_threads.to_string()])
        .args(["-k", &args.mmseqs_args.k.to_string()])
        .args(["--k-score", &args.mmseqs_args.k_score.to_string()])
        .args([
            "--min-ungapped-score",
            &args.mmseqs_args.min_ungapped_score.to_string(),
        ])
        .args(["--max-seqs", &args.mmseqs_args.max_seqs.to_string()])
        .args(["-e", &effective_e_value.to_string()])
        // the '-a' argument enables alignment backtraces in mmseqs2
        // it is required to get start positions for alignments
        .args(["-a", "1"])
        .run()?;

    Command::new("mmseqs")
        .arg("convertalis")
        .arg(&args.prep_dir.mmseqs_query_db_path())
        .arg(&args.prep_dir.mmseqs_target_db_path())
        .arg(&args.prep_dir.mmseqs_align_db_path())
        .arg(&args.prep_dir.mmseqs_align_tsv_path())
        .args(["--threads", &args.common_args.num_threads.to_string()])
        .args([
            "--format-output",
            "qheader,theader,qstart,qend,tstart,tend,evalue",
        ])
        .run()?;
    Ok(())
}

pub fn seed(args: &SeedArgs) -> anyhow::Result<SeedMap> {
    run_mmseqs_search(args)?;

    let profile_seeds_by_name =
        build_alignment_seeds(args).context("failed to build alignment seeds")?;

    let mut seeds_out = args
        .seeds_path
        .open(true)
        .context("failed to create alignment seeds file")?;

    write!(
        seeds_out,
        "{}",
        serde_json::to_string(&profile_seeds_by_name)?
    )
    .context("failed to write alignment seeds")?;

    Ok(profile_seeds_by_name)
}

pub fn map_p7_to_mmseqs_profiles(
    p7_profiles: &[Profile],
    args: &SeedArgs,
) -> anyhow::Result<HashMap<String, Vec<usize>>> {
    let mmseqs_consensus_map = extract_mmseqs_profile_consensus_sequences(args)?;

    let mut mmseqs_to_p7_idx_by_name: HashMap<String, Vec<usize>> = HashMap::new();

    for p7_profile in p7_profiles.iter() {
        let name = &p7_profile.name;
        let mmseqs_consensus = mmseqs_consensus_map.get(name).unwrap();
        let p7_consensus = Sequence::from_utf8(&p7_profile.consensus_sequence_bytes_utf8[1..])?;
        let trace = needleman_wunsch(mmseqs_consensus, &p7_consensus);

        let mut mmseqs_to_p7: Vec<usize> = vec![0; mmseqs_consensus.length + 1];

        let mut mmseqs_idx: usize = 0;
        let mut p7_idx: usize = 0;
        for step in &trace {
            match step {
                SimpleTraceStep::Diagonal => {
                    mmseqs_idx += 1;
                    p7_idx += 1;
                }
                SimpleTraceStep::Up => {
                    mmseqs_idx += 1;
                }
                SimpleTraceStep::Left => {
                    p7_idx += 1;
                }
            }
            mmseqs_to_p7[mmseqs_idx] = p7_idx;
        }

        // this debug assert should guarantee that the NW
        // alignment fully covered both consensus sequences
        debug_assert_eq!(mmseqs_idx, mmseqs_consensus.length);
        debug_assert_eq!(p7_idx, p7_consensus.length);

        mmseqs_to_p7_idx_by_name.insert(name.clone(), mmseqs_to_p7);
    }

    Ok(mmseqs_to_p7_idx_by_name)
}

pub fn extract_mmseqs_profile_consensus_sequences(
    args: &SeedArgs,
) -> anyhow::Result<HashMap<String, Sequence>> {
    let mut offsets_and_lengths: Vec<(usize, usize)> = vec![];
    let mut names: Vec<String> = vec![];

    // query_db_h_index -> offsets & lengths of query headers
    // query_db_h       -> query headers
    // query_db_index   -> offsets & lengths of binary profiles
    // query_db         -> binary profiles (lines of 23 or 25, depending on mmseqs version)

    let query_db_h_index_file = File::open(args.prep_dir.mmseqs_query_db_h_index_path())
        .context("failed to open queryDB_h.index")?;

    let reader = BufReader::new(query_db_h_index_file);
    for line in reader.lines() {
        match line {
            Ok(line) => {
                let tokens: Vec<&str> = line.split_whitespace().collect();

                let offset = tokens[1].parse::<usize>()?;
                let length = tokens[2].parse::<usize>()?;
                offsets_and_lengths.push((offset, length));
            }
            Err(e) => {
                return Err(e).context("failed to parse line in queryDB_h.index");
            }
        }
    }

    let mut query_db_h_file =
        File::open(args.prep_dir.mmseqs_query_db_h_path()).context("failed to open queryDB_h")?;

    for (offset, length) in &offsets_and_lengths {
        let mut buffer = vec![0; *length];
        query_db_h_file.seek(SeekFrom::Start(*offset as u64))?;
        query_db_h_file.read_exact(&mut buffer)?;

        let mut name: Option<String> = None;
        for (buf_idx, byte) in buffer.iter().enumerate() {
            if byte.is_ascii_whitespace() {
                name = Some(
                    std::str::from_utf8(&buffer[0..buf_idx])
                        .context("failed to create profile name string")?
                        .to_string(),
                );
                break;
            }
        }

        match name {
            Some(name) => names.push(name),
            None => {
                panic!()
            }
        }
    }

    let query_db_index_file = File::open(args.prep_dir.mmseqs_query_db_index_path())
        .context("failed to open queryDB.index")?;

    let reader = BufReader::new(query_db_index_file);
    for line in reader.lines() {
        match line {
            Ok(line) => {
                let tokens: Vec<&str> = line.split_whitespace().collect();

                let line_idx = tokens[0].parse::<usize>()?;
                let offset = tokens[1].parse::<usize>()?;
                let length = tokens[2].parse::<usize>()?;
                offsets_and_lengths[line_idx] = (offset, length);
            }
            Err(e) => {
                return Err(e).context("failed to parse line in queryDB.index");
            }
        }
    }

    let mut sequence_map: HashMap<String, Sequence> = HashMap::new();

    let mut query_db_file =
        File::open(args.prep_dir.mmseqs_query_db_path()).context("failed to open queryDB")?;

    for (seq_idx, (offset, length)) in offsets_and_lengths.iter().enumerate() {
        let mut buffer = vec![0; *length];
        query_db_file.seek(SeekFrom::Start(*offset as u64))?;
        query_db_file.read_exact(&mut buffer)?;

        let mut consensus_digital_bytes: Vec<u8> = vec![];

        for byte_chunk in buffer.chunks(25) {
            if byte_chunk.len() == 25 {
                consensus_digital_bytes.push(byte_chunk[21]);
            }
        }

        sequence_map.insert(
            names[seq_idx].clone(),
            Sequence::from_digital(&consensus_digital_bytes)?,
        );
    }

    Ok(sequence_map)
}

#[derive(Error, Debug)]
#[error("no profile to profile map for: {profile_name}")]
pub struct ProfilesNotMappedError {
    pub profile_name: String,
}

// TODO: move this elsewhere
#[derive(Error, Debug)]
#[error("invalid file format: {format}")]
pub struct InvalidFileFormatError {
    pub format: FileFormat,
}

pub fn build_alignment_seeds(args: &SeedArgs) -> anyhow::Result<SeedMap> {
    let mut seed_map: SeedMap = HashMap::new();

    let mmseqs_align_file = File::open(args.prep_dir.mmseqs_align_tsv_path()).context(format!(
        "couldn't open mmseqs align file at: {}",
        &args.prep_dir.mmseqs_align_tsv_path().to_string_lossy()
    ))?;

    let align_reader = BufReader::new(mmseqs_align_file);

    for line in align_reader.lines().flatten() {
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
