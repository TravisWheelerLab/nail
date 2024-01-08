use std::collections::HashMap;
use std::fs::File;
use std::io::Read;
use std::io::{stdout, BufWriter, Write};
use std::path::PathBuf;
use std::sync::Mutex;

use anyhow::Context;
use clap::Args;
use rayon::iter::{IntoParallelIterator, ParallelIterator};
use thiserror::Error;

use crate::args::{guess_query_format_from_query_file, FileFormat};
use crate::cli::CommonArgs;
use crate::extension_traits::PathBufExt;
use crate::pipeline::prep::{build_hmm_from_fasta, build_hmm_from_stockholm};
use crate::pipeline::seed::SeedMap;

use libnail::align::structs::{
    Alignment, CloudBoundGroup, CloudMatrixLinear, CloudSearchParams, DpMatrixSparse, RowBounds,
    ScoreParams, Seed, Trace,
};
use libnail::align::{
    backward, cloud_search_backward, cloud_search_forward, forward, null1_score,
    null2_score_bounded, optimal_accuracy, posterior, traceback,
};
use libnail::structs::hmm::parse_hmms_from_p7hmm_file;
use libnail::structs::{Profile, Sequence};

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
    pub nail_args: NailArgs,

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

    align_loop(args, profiles, targets, seed_map)?;

    Ok(())
}

pub fn align_loop(
    args: &AlignArgs,
    mut profiles: Vec<Profile>,
    targets: Vec<Sequence>,
    seed_map: SeedMap,
) -> anyhow::Result<()> {
    let tab_results_writer: Mutex<BufWriter<File>> =
        Mutex::new(args.output_args.tsv_results_path.open(true)?);
    let ali_results_writer: Mutex<Box<dyn Write + Send>> = match args.output_args.ali_results_path {
        Some(ref path) => Mutex::new(Box::new(path.open(true)?)),
        None => Mutex::new(Box::new(stdout())),
    };

    let dp = AlignmentStructs::default();

    let num_targets = match args.nail_args.target_database_size {
        Some(size) => size,
        None => targets.len(),
    };

    let score_params = ScoreParams::new(num_targets);

    let cloud_search_params = CloudSearchParams {
        gamma: args.nail_args.gamma,
        alpha: args.nail_args.alpha,
        beta: args.nail_args.beta,
    };

    let mut target_map: HashMap<String, Sequence> = HashMap::new();
    for target in targets {
        target_map.insert(target.name.clone(), target);
    }

    let mut profile_seeds_pairs: Vec<(&mut Profile, &Vec<Seed>)> = vec![];

    for profile in profiles.iter_mut() {
        match seed_map.get(&profile.name) {
            Some(seeds) => profile_seeds_pairs.push((profile, seeds)),
            None => {
                continue;
            }
        }
    }

    #[derive(Default, Clone)]
    struct AlignmentStructs {
        cloud_matrix: CloudMatrixLinear,
        forward_bounds: CloudBoundGroup,
        backward_bounds: CloudBoundGroup,
        forward_matrix: DpMatrixSparse,
        backward_matrix: DpMatrixSparse,
        posterior_matrix: DpMatrixSparse,
        optimal_matrix: DpMatrixSparse,
    }

    profile_seeds_pairs.into_par_iter().for_each_with(
        (dp, score_params),
        |(dp, score_params), (profile, seeds)| {
            for seed in seeds {
                let target = target_map.get(&seed.target_name).unwrap();
                profile.configure_for_target_length(target.length);

                dp.cloud_matrix.reuse(profile.length);
                dp.forward_bounds.reuse(target.length, profile.length);
                dp.backward_bounds.reuse(target.length, profile.length);

                cloud_search_forward(
                    profile,
                    target,
                    seed,
                    &mut dp.cloud_matrix,
                    &cloud_search_params,
                    &mut dp.forward_bounds,
                );

                cloud_search_backward(
                    profile,
                    target,
                    seed,
                    &mut dp.cloud_matrix,
                    &cloud_search_params,
                    &mut dp.backward_bounds,
                );

                if dp.forward_bounds.max_anti_diagonal_idx
                    < dp.backward_bounds.min_anti_diagonal_idx
                {
                    dp.forward_bounds.fill_rectangle(
                        seed.target_start,
                        seed.profile_start,
                        seed.target_end,
                        seed.profile_end,
                    );
                } else {
                    CloudBoundGroup::join_merge(&mut dp.forward_bounds, &dp.backward_bounds);
                }

                if !dp.forward_bounds.valid() {
                    dp.forward_bounds.fill_rectangle(
                        seed.target_start,
                        seed.profile_start,
                        seed.target_end,
                        seed.profile_end,
                    );
                }

                dp.forward_bounds.trim_wings();

                let mut row_bounds = RowBounds::new(&dp.forward_bounds);

                if !row_bounds.valid() {
                    dp.forward_bounds.fill_rectangle(
                        seed.target_start,
                        seed.profile_start,
                        seed.target_end,
                        seed.profile_end,
                    );
                    row_bounds = RowBounds::new(&dp.forward_bounds);
                }

                dp.forward_matrix
                    .reuse(target.length, profile.length, &row_bounds);
                dp.backward_matrix
                    .reuse(target.length, profile.length, &row_bounds);
                dp.posterior_matrix
                    .reuse(target.length, profile.length, &row_bounds);
                dp.optimal_matrix
                    .reuse(target.length, profile.length, &row_bounds);

                // we use the forward score to compute the final bit score (later)
                score_params.forward_score_nats =
                    forward(profile, target, &mut dp.forward_matrix, &row_bounds);

                backward(profile, target, &mut dp.backward_matrix, &row_bounds);

                posterior(
                    profile,
                    &dp.forward_matrix,
                    &dp.backward_matrix,
                    &mut dp.posterior_matrix,
                    &row_bounds,
                );

                score_params.null_score_nats = null1_score(target.length);
                score_params.bias_correction_score_nats =
                    null2_score_bounded(&dp.posterior_matrix, profile, target, &row_bounds);

                optimal_accuracy(
                    profile,
                    &dp.posterior_matrix,
                    &mut dp.optimal_matrix,
                    &row_bounds,
                );

                let mut trace = Trace::new(target.length, profile.length);
                traceback(
                    profile,
                    &dp.posterior_matrix,
                    &dp.optimal_matrix,
                    &mut trace,
                    row_bounds.target_end,
                );

                let alignment = Alignment::from_trace(&trace, profile, target, score_params);

                if alignment.evalue <= args.output_args.evalue_threshold {
                    match tab_results_writer.lock() {
                        Ok(mut writer) => writeln!(writer, "{}", alignment.tab_string())
                            .expect("failed to write tabular output"),
                        Err(_) => panic!("tabular results writer mutex was poisoned"),
                    };

                    match ali_results_writer.lock() {
                        Ok(mut writer) => {
                            writeln!(writer, "{}", alignment.ali_string())
                                .expect("failed to write alignment output");
                        }
                        Err(_) => panic!("alignment results writer mutex was poisoned"),
                    }
                }
            }
        },
    );

    Ok(())
}
