mod serial_bounded;
pub use serial_bounded::*;
mod serial_full;
pub use serial_full::*;
mod threaded_bounded;
pub use threaded_bounded::*;
mod threaded_full;
pub use threaded_full::*;

use std::fs::File;
use std::io::Read;
use std::path::PathBuf;

use crate::args::{guess_query_format_from_query_file, FileFormat};
use crate::pipeline::prep::{build_hmm_from_fasta, build_hmm_from_stockholm};
use crate::pipeline::seed::SeedMap;

use libnail::structs::hmm::parse_hmms_from_p7hmm_file;

use crate::cli::CommonArgs;
use anyhow::Context;
use clap::Args;
use libnail::structs::{Profile, Sequence};
use thiserror::Error;

#[derive(Error, Debug)]
#[error("no profile with name: {profile_name}")]
pub struct ProfileNotFoundError {
    profile_name: String,
}

#[derive(Error, Debug)]
#[error("no target with name: {target_name}")]
pub struct TargetNotFoundError {
    target_name: String,
}

#[derive(Args, Debug, Clone)]
pub struct NailArgs {
    /// Override the target database size (number of sequences) used for E-value calculation
    #[arg(short = 'Z', value_name = "N")]
    pub target_database_size: Option<usize>,
    /// Pruning parameter alpha
    #[arg(short = 'A', default_value_t = 12.0, value_name = "F")]
    pub alpha: f32,
    /// Pruning parameter beta
    #[arg(short = 'B', default_value_t = 20.0, value_name = "F")]
    pub beta: f32,
    /// Pruning parameter gamma
    #[arg(short = 'G', default_value_t = 5, value_name = "N")]
    pub gamma: usize,
    /// Compute the full dynamic programming matrices during alignment
    #[arg(long, action)]
    pub full_dp: bool,
}

#[derive(Args, Debug, Clone)]
pub struct AlignOutputArgs {
    /// Only report hits with an E-value below this value
    #[arg(short = 'E', default_value_t = 10.0, value_name = "F")]
    pub evalue_threshold: f64,
    /// Where to place tabular output
    #[arg(
        short = 'T',
        long = "tab-output",
        default_value = "results.tsv",
        value_name = "path"
    )]
    pub tsv_results_path: PathBuf,
    /// Where to place alignment output
    #[arg(short = 'O', long = "output", value_name = "path")]
    pub ali_results_path: Option<PathBuf>,
}

#[derive(Debug, Args)]
pub struct AlignArgs {
    /// Query file
    #[arg(value_name = "QUERY.[fasta:hmm:sto]")]
    pub query_path: PathBuf,
    /// Target file
    #[arg(value_name = "TARGET.fasta")]
    pub target_path: PathBuf,
    /// Alignment seeds from running nail seed (or elsewhere)
    #[arg(value_name = "SEEDS.json")]
    pub seeds_path: PathBuf,

    /// Arguments that are passed to nail functions
    #[command(flatten)]
    pub libnail_args: NailArgs,

    /// Arguments that control output options
    #[command(flatten)]
    pub output_args: AlignOutputArgs,

    /// Arguments that are common across all nail subcommands
    #[command(flatten)]
    pub common_args: CommonArgs,
}

pub fn align(
    args: &AlignArgs,
    profiles: Option<Vec<Profile>>,
    seed_map: Option<SeedMap>,
) -> anyhow::Result<()> {
    let profiles = match profiles {
        // if we happened to run the seed step before
        // this, the profiles will be passed in
        Some(profiles) => profiles,
        None => {
            let query_format = guess_query_format_from_query_file(&args.query_path)?;
            let hmm_path = match query_format {
                FileFormat::Fasta => {
                    let hmm_path = args.query_path.with_extension("hmm");
                    build_hmm_from_fasta(
                        &args.query_path,
                        &hmm_path,
                        args.common_args.num_threads,
                    )?;
                    hmm_path
                }
                FileFormat::Stockholm => {
                    let hmm_path = args.query_path.with_extension("hmm");
                    build_hmm_from_stockholm(
                        &args.query_path,
                        &hmm_path,
                        args.common_args.num_threads,
                    )?;
                    hmm_path
                }
                FileFormat::Hmm => args.query_path.clone(),
                FileFormat::Unset => {
                    // TODO: real error
                    panic!("query format is unset in call to align()");
                }
            };

            let hmms = parse_hmms_from_p7hmm_file(hmm_path)?;

            hmms.iter().map(Profile::new).collect()
        }
    };

    let seed_map = match seed_map {
        // if we happened to run the seed step before
        // this, the seeds will be passed in
        Some(seed_map) => seed_map,
        None => {
            let mut seeds_string = String::new();
            File::open(&args.seeds_path)
                .context(format!(
                    "failed to open alignment seeds file: {}",
                    &args.seeds_path.to_string_lossy(),
                ))?
                .read_to_string(&mut seeds_string)
                .context(format!(
                    "failed to read alignment seeds file: {}",
                    &args.seeds_path.to_string_lossy(),
                ))?;

            serde_json::from_str(&seeds_string).context(format!(
                "failed to parse alignment seeds file: {}",
                &args.seeds_path.to_string_lossy(),
            ))?
        }
    };

    let targets = Sequence::amino_from_fasta(&args.target_path)?;

    // this is how we tell rayon how many threads to use
    rayon::ThreadPoolBuilder::new()
        .num_threads(args.common_args.num_threads)
        .build_global()
        .unwrap();

    if args.libnail_args.full_dp {
        align_threaded_full(args, profiles, targets, seed_map)?;
    } else {
        align_threaded_bounded(args, profiles, targets, seed_map)?;
    }
    // }

    Ok(())
}
