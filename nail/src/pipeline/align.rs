use std::collections::HashMap;
use std::fs::File;
use std::io::Read;
use std::io::{stdout, BufWriter, Write};
use std::path::PathBuf;
use std::sync::{Arc, Mutex};

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
    Alignment, AntiDiagonalBounds, CloudMatrixLinear, CloudSearchParams, DpMatrixSparse, RowBounds,
    ScoreParams, Seed, Trace,
};
use libnail::align::{
    backward, cloud_search_backward, cloud_search_forward, composition_bias_score, forward,
    length_bias_score, optimal_accuracy, posterior, traceback,
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
    /// The P-value threshold for promoting hits past cloud search
    #[arg(long = "cloud-thresh", default_value_t = 1e-3, value_name = "F")]
    pub cloud_pvalue_threshold: f64,
    /// The P-value threshold for promoting hits past forward
    #[arg(long = "forward-thresh", default_value_t = 1e-4, value_name = "F")]
    pub forward_pvalue_threshold: f64,
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

#[derive(Clone)]
struct OutputData {
    tab_results_writer: Arc<Mutex<BufWriter<File>>>,
    ali_results_writer: Arc<Mutex<Box<dyn Write + Send>>>,
}

#[derive(Default, Clone)]
struct CloudSearchData {
    cloud_search_params: CloudSearchParams,
    cloud_matrix: CloudMatrixLinear,
    forward_bounds: AntiDiagonalBounds,
    backward_bounds: AntiDiagonalBounds,
    row_bounds: RowBounds,
}

#[derive(Default, Clone)]
struct AlignmentData {
    forward_matrix: DpMatrixSparse,
    backward_matrix: DpMatrixSparse,
    posterior_matrix: DpMatrixSparse,
    optimal_matrix: DpMatrixSparse,
}

#[derive(Clone)]
struct ThreadData {
    output: OutputData,
    cloud_search: CloudSearchData,
    alignment: AlignmentData,
    target_map: Arc<HashMap<String, Sequence>>,
    score_params: ScoreParams,
    evalue_threshold: f64,
    cloud_pvalue_threshold: f64,
    forward_pvalue_threshold: f64,
    full_dp: bool,
}

impl ThreadData {
    pub fn new(args: &AlignArgs, targets: Vec<Sequence>) -> anyhow::Result<Self> {
        let tab_writer = Mutex::new(args.output_args.tsv_results_path.open(true)?);

        // TODO: is there a better way to do this without a trait object?
        // note: I think the type is required here to tell
        //       the compiler that this is a trait object
        let ali_writer: Mutex<Box<dyn Write + Send>> = match args.output_args.ali_results_path {
            Some(ref path) => Mutex::new(Box::new(path.open(true)?)),
            None => Mutex::new(Box::new(stdout())),
        };

        let output = OutputData {
            tab_results_writer: Arc::new(tab_writer),
            ali_results_writer: Arc::new(ali_writer),
        };

        let num_targets = match args.nail_args.target_database_size {
            Some(size) => size,
            None => targets.len(),
        };

        let score_params = ScoreParams::new(num_targets);

        let mut target_map: HashMap<String, Sequence> = HashMap::new();
        for target in targets {
            target_map.insert(target.name.clone(), target);
        }

        let target_map = Arc::new(target_map);

        Ok(Self {
            output,
            cloud_search: CloudSearchData {
                cloud_search_params: CloudSearchParams {
                    gamma: args.nail_args.gamma,
                    alpha: args.nail_args.alpha,
                    beta: args.nail_args.beta,
                },
                ..Default::default()
            },
            alignment: AlignmentData::default(),
            target_map,
            score_params,
            evalue_threshold: args.output_args.evalue_threshold,
            cloud_pvalue_threshold: args.nail_args.cloud_pvalue_threshold,
            forward_pvalue_threshold: args.nail_args.forward_pvalue_threshold,
            full_dp: args.nail_args.full_dp,
        })
    }
}

struct ProfileSeedsPair<'a> {
    profile: &'a mut Profile,
    seeds: &'a Vec<Seed>,
}

pub fn align(
    args: &AlignArgs,
    profiles: Option<Vec<Profile>>,
    seed_map: Option<SeedMap>,
) -> anyhow::Result<()> {
    let mut profiles = match profiles {
        // if we happened to run the seed step before
        // this, the profiles will be passed in
        Some(profiles) => profiles,
        None => {
            let query_format = guess_query_format_from_query_file(&args.query_path)?;
            let hmm_path = match query_format {
                FileFormat::Fasta => {
                    let hmm_path = args.query_path.with_extension("hmm");
                    let query_dir = args.query_path.parent().with_context(|| {
                        format!(
                            "failed to resolve query directory from path: {}",
                            args.query_path.to_string_lossy()
                        )
                    })?;
                    let temp_hmm_path = query_dir.join("tmp.hmm");
                    let temp_fasta_path = query_dir.join("tmp.fa");

                    build_hmm_from_fasta(
                        &args.query_path,
                        &temp_fasta_path,
                        &hmm_path,
                        &temp_hmm_path,
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

    // the thread data is everything that is cloned per thread
    let thread_data = ThreadData::new(args, targets)?;

    match thread_data.output.tab_results_writer.lock() {
        Ok(mut writer) => {
            writeln!(writer, "{}", Alignment::TAB_HEADER)
        }
        Err(_) => panic!("tabular results writer mutex was poisoned"),
    }?;

    let mut profile_seeds_pairs: Vec<ProfileSeedsPair> = vec![];

    for profile in profiles.iter_mut() {
        match seed_map.get(&profile.name) {
            Some(seeds) => profile_seeds_pairs.push(ProfileSeedsPair { profile, seeds }),
            None => {
                continue;
            }
        }
    }

    // this is how we tell rayon how many threads to use
    rayon::ThreadPoolBuilder::new()
        .num_threads(args.common_args.num_threads)
        .build_global()
        .unwrap();

    profile_seeds_pairs
        .into_par_iter()
        .for_each_with(thread_data, align_seeds);

    Ok(())
}

fn cloud_search(
    profile: &Profile,
    target: &Sequence,
    seed: &Seed,
    dp: &mut CloudSearchData,
) -> f32 {
    dp.cloud_matrix.reuse(profile.length);
    dp.forward_bounds.reuse(target.length, profile.length);
    dp.backward_bounds.reuse(target.length, profile.length);
    dp.row_bounds.reuse(target.length);

    let forward_scores = cloud_search_forward(
        profile,
        target,
        seed,
        &mut dp.cloud_matrix,
        &dp.cloud_search_params,
        &mut dp.forward_bounds,
    );

    let backward_scores = cloud_search_backward(
        profile,
        target,
        seed,
        &mut dp.cloud_matrix,
        &dp.cloud_search_params,
        &mut dp.backward_bounds,
    );

    // this approximates the score for the forward
    // cloud that extends past the seed end point
    let disjoint_forward_score = forward_scores.max_score - forward_scores.max_score_within;

    // this approximates the score for the backward
    // cloud that extends past the seed start point
    let disjoint_backward_score = backward_scores.max_score - backward_scores.max_score_within;

    // this approximates the score of the intersection
    // of the forward and backward clouds
    let intersection_score = forward_scores
        .max_score_within
        .max(backward_scores.max_score_within);

    let cloud_score_nats = intersection_score + disjoint_forward_score + disjoint_backward_score;
    let cloud_score_bits = cloud_score_nats / std::f32::consts::LN_2;

    let bounds_intersect =
        dp.forward_bounds.max_anti_diagonal_idx >= dp.backward_bounds.min_anti_diagonal_idx;

    if bounds_intersect {
        AntiDiagonalBounds::join_merge(&mut dp.forward_bounds, &dp.backward_bounds);
        if !dp.forward_bounds.valid() {
            dp.forward_bounds.fill_rectangle(
                seed.target_start,
                seed.profile_start,
                seed.target_end,
                seed.profile_end,
            );
        }

        dp.forward_bounds.trim_wings();

        dp.row_bounds
            .fill_from_anti_diagonal_bounds(&dp.forward_bounds);

        if !dp.row_bounds.valid() {
            dp.row_bounds.fill_rectangle(
                seed.target_start,
                seed.profile_start,
                seed.target_end,
                seed.profile_end,
            );
        }
    } else {
        dp.row_bounds.fill_rectangle(
            seed.target_start,
            seed.profile_start,
            seed.target_end,
            seed.profile_end,
        );
    }

    cloud_score_bits
}

fn align_seeds(data: &mut ThreadData, pair: ProfileSeedsPair) {
    let alignment_data = &mut data.alignment;
    let cloud_search_data = &mut data.cloud_search;
    let profile_length = pair.profile.length;

    for seed in pair.seeds {
        let target = data.target_map.get(&seed.target_name).unwrap();
        pair.profile.configure_for_target_length(target.length);

        let cloud_score_bits = if data.full_dp {
            cloud_search_data
                .row_bounds
                .fill_rectangle(1, 1, target.length, profile_length);
            f32::MAX
        } else {
            cloud_search(pair.profile, target, seed, cloud_search_data)
        };

        let cloud_pvalue = (-pair.profile.forward_lambda as f64
            * (cloud_score_bits as f64 - pair.profile.forward_tau as f64))
            .exp();

        if cloud_pvalue >= data.cloud_pvalue_threshold {
            continue;
        }

        alignment_data.forward_matrix.reuse(
            target.length,
            profile_length,
            &cloud_search_data.row_bounds,
        );
        alignment_data.backward_matrix.reuse(
            target.length,
            profile_length,
            &cloud_search_data.row_bounds,
        );
        alignment_data.posterior_matrix.reuse(
            target.length,
            profile_length,
            &cloud_search_data.row_bounds,
        );
        alignment_data.optimal_matrix.reuse(
            target.length,
            profile_length,
            &cloud_search_data.row_bounds,
        );

        // we use the forward score to compute the final bit score (later)
        data.score_params.forward_score_nats = forward(
            pair.profile,
            target,
            &mut alignment_data.forward_matrix,
            &cloud_search_data.row_bounds,
        );

        let forward_pvalue = (-pair.profile.forward_lambda as f64
            * ((data.score_params.forward_score_nats / std::f32::consts::LN_2) as f64
                - pair.profile.forward_tau as f64))
            .exp();

        if forward_pvalue >= data.forward_pvalue_threshold {
            continue;
        }

        backward(
            pair.profile,
            target,
            &mut alignment_data.backward_matrix,
            &cloud_search_data.row_bounds,
        );

        posterior(
            pair.profile,
            &alignment_data.forward_matrix,
            &alignment_data.backward_matrix,
            &mut alignment_data.posterior_matrix,
            &cloud_search_data.row_bounds,
        );

        data.score_params.length_bias_score_nats = length_bias_score(target.length);

        data.score_params.composition_bias_score_nats = composition_bias_score(
            &alignment_data.posterior_matrix,
            pair.profile,
            target,
            &cloud_search_data.row_bounds,
        );

        optimal_accuracy(
            pair.profile,
            &alignment_data.posterior_matrix,
            &mut alignment_data.optimal_matrix,
            &cloud_search_data.row_bounds,
        );

        let mut trace = Trace::new(target.length, profile_length);
        traceback(
            pair.profile,
            &alignment_data.posterior_matrix,
            &alignment_data.optimal_matrix,
            &mut trace,
            cloud_search_data.row_bounds.target_end,
        );

        let mut alignment = Alignment::from_trace(&trace, pair.profile, target, &data.score_params);

        alignment.cell_fraction = Some(
            cloud_search_data.row_bounds.num_cells() as f32
                / (target.length * profile_length) as f32,
        );

        if alignment.evalue <= data.evalue_threshold {
            match data.output.tab_results_writer.lock() {
                Ok(mut writer) => writeln!(writer, "{}", alignment.tab_string())
                    .expect("failed to write tabular output"),
                Err(_) => panic!("tabular results writer mutex was poisoned"),
            };

            match data.output.ali_results_writer.lock() {
                Ok(mut writer) => {
                    writeln!(writer, "{}", alignment.ali_string())
                        .expect("failed to write alignment output");
                }
                Err(_) => panic!("alignment results writer mutex was poisoned"),
            }
        }
    }
}
