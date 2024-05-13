use std::collections::HashMap;
use std::io::{stdout, Write};
use std::path::PathBuf;
use std::sync::{Arc, Mutex};

use clap::Args;
use libnail::output::output_tabular::{Field, TableFormat};
use rayon::iter::{IntoParallelIterator, ParallelIterator};
use thiserror::Error;

use crate::cli::CommonArgs;
use crate::database::Database;
use crate::extension_traits::PathBufExt;

use libnail::align::structs::{
    Alignment, AlignmentBuilder, AntiDiagonalBounds, CloudMatrixLinear, DpMatrixSparse, RowBounds,
    Seed, Trace,
};
use libnail::align::{
    backward, cloud_score, cloud_search_backward, cloud_search_forward, forward, null_one_score,
    null_two_score, optimal_accuracy, p_value, posterior, traceback, CloudSearchParams,
};
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

pub const DEFAULT_COLUMNS: [Field; 10] = [
    Field::Target,
    Field::Query,
    Field::TargetStart,
    Field::TargetEnd,
    Field::QueryStart,
    Field::QueryEnd,
    Field::Score,
    Field::CompBias,
    Field::Evalue,
    Field::CellFrac,
];

pub trait SeedStep: Clone {
    fn run(&self, profile: &Profile, target: &Sequence) -> Option<Seed>;
}

pub trait CloudSearchStep: Clone {
    fn run(&mut self, profile: &Profile, target: &Sequence, seed: &Seed) -> Option<&RowBounds>;
}

pub trait AlignStep: Clone {
    fn run(
        &mut self,
        profile: &mut Profile,
        target: &Sequence,
        bounds: &RowBounds,
    ) -> Option<Alignment>;
}

#[derive(Default, Clone)]
pub struct DefaultSeedStep {
    seeds: HashMap<String, HashMap<String, Seed>>,
}

impl DefaultSeedStep {
    pub fn new(seeds: HashMap<String, HashMap<String, Seed>>) -> Self {
        DefaultSeedStep { seeds }
    }
}

impl SeedStep for DefaultSeedStep {
    fn run(&self, profile: &Profile, target: &Sequence) -> Option<Seed> {
        Some(self.seeds.get(&profile.name)?.get(&target.name)?.clone())
    }
}

#[derive(Default, Clone)]
pub struct DefaultCloudSearchStep {
    cloud_matrix: CloudMatrixLinear,
    forward_bounds: AntiDiagonalBounds,
    reverse_bounds: AntiDiagonalBounds,
    row_bounds: RowBounds,
    params: CloudSearchParams,
    p_value_threshold: f64,
}

impl DefaultCloudSearchStep {
    pub fn new(args: &AlignArgs) -> Self {
        Self {
            params: CloudSearchParams {
                gamma: args.nail_args.gamma,
                alpha: args.nail_args.alpha,
                beta: args.nail_args.beta,
            },
            p_value_threshold: args.nail_args.cloud_pvalue_threshold,
            ..Default::default()
        }
    }
}

impl CloudSearchStep for DefaultCloudSearchStep {
    fn run(&mut self, profile: &Profile, target: &Sequence, seed: &Seed) -> Option<&RowBounds> {
        self.cloud_matrix.reuse(profile.length);
        self.forward_bounds.reuse(target.length, profile.length);
        self.reverse_bounds.reuse(target.length, profile.length);
        self.row_bounds.reuse(target.length);

        let forward_scores = cloud_search_forward(
            profile,
            target,
            seed,
            &mut self.cloud_matrix,
            &self.params,
            &mut self.forward_bounds,
        );

        let reverse_scores = cloud_search_backward(
            profile,
            target,
            seed,
            &mut self.cloud_matrix,
            &self.params,
            &mut self.reverse_bounds,
        );

        let cloud_score = cloud_score(&forward_scores, &reverse_scores);

        let cloud_p_value = p_value(cloud_score, profile.forward_lambda, profile.forward_tau);

        if cloud_p_value >= self.p_value_threshold {
            return None;
        }

        let bounds_intersect =
            self.forward_bounds.max_anti_diagonal_idx >= self.reverse_bounds.min_anti_diagonal_idx;

        // TODO: clean up this mess
        if bounds_intersect {
            AntiDiagonalBounds::join_merge(&mut self.forward_bounds, &self.reverse_bounds);
            if !self.forward_bounds.valid() {
                self.forward_bounds.fill_rectangle(
                    seed.target_start,
                    seed.profile_start,
                    seed.target_end,
                    seed.profile_end,
                );
            }

            self.forward_bounds.square_corners();
            self.forward_bounds.trim_wings();

            self.row_bounds
                .fill_from_anti_diagonal_bounds(&self.forward_bounds);

            if !self.row_bounds.valid() {
                self.row_bounds.fill_rectangle(
                    seed.target_start,
                    seed.profile_start,
                    seed.target_end,
                    seed.profile_end,
                );
            }
        } else {
            self.row_bounds.fill_rectangle(
                seed.target_start,
                seed.profile_start,
                seed.target_end,
                seed.profile_end,
            );
        }

        Some(&self.row_bounds)
    }
}

#[derive(Default, Clone)]
pub struct DefaultAlignStep {
    forward_matrix: DpMatrixSparse,
    backward_matrix: DpMatrixSparse,
    posterior_matrix: DpMatrixSparse,
    optimal_matrix: DpMatrixSparse,
    forward_pvalue_threshold: f64,
    target_count: usize,
    e_value_threshold: f64,
}

impl DefaultAlignStep {
    pub fn new(args: &AlignArgs, target_count: usize) -> Self {
        Self {
            target_count: match args.nail_args.target_database_size {
                Some(size) => size,
                None => target_count,
            },
            forward_pvalue_threshold: args.nail_args.forward_pvalue_threshold,
            e_value_threshold: args.output_args.evalue_threshold,
            ..Default::default()
        }
    }
}

impl AlignStep for DefaultAlignStep {
    fn run(
        &mut self,
        profile: &mut Profile,
        target: &Sequence,
        bounds: &RowBounds,
    ) -> Option<Alignment> {
        // configuring for the target length
        // adjusts special state transitions
        profile.configure_for_target_length(target.length);

        self.forward_matrix
            .reuse(target.length, profile.length, bounds);
        self.backward_matrix
            .reuse(target.length, profile.length, bounds);
        self.posterior_matrix
            .reuse(target.length, profile.length, bounds);
        self.optimal_matrix
            .reuse(target.length, profile.length, bounds);

        // we use the forward score to compute the final bit score (later)
        let forward_score = forward(profile, target, &mut self.forward_matrix, bounds).to_bits();

        // for now we compute the P-value for filtering purposes
        let forward_pvalue = p_value(forward_score, profile.forward_lambda, profile.forward_tau);

        if forward_pvalue >= self.forward_pvalue_threshold {
            return None;
        }

        backward(profile, target, &mut self.backward_matrix, bounds);

        posterior(
            profile,
            &self.forward_matrix,
            &self.backward_matrix,
            &mut self.posterior_matrix,
            bounds,
        );

        optimal_accuracy(
            profile,
            &self.posterior_matrix,
            &mut self.optimal_matrix,
            bounds,
        );

        let mut trace = Trace::new(target.length, profile.length);
        traceback(
            profile,
            &self.posterior_matrix,
            &self.optimal_matrix,
            &mut trace,
            bounds.target_end,
        );

        let null_one = null_one_score(target.length);
        let null_two = null_two_score(&self.posterior_matrix, profile, target, bounds);

        let cell_fraction = bounds.num_cells() as f32 / (profile.length * target.length) as f32;

        let alignment = AlignmentBuilder::new(&trace)
            .with_profile(profile)
            .with_target(target)
            .with_target_count(self.target_count)
            .with_forward_score(forward_score)
            .with_null_one(null_one)
            .with_null_two(null_two)
            .with_cell_fraction(cell_fraction)
            .build()
            .unwrap();

        match alignment.e_value {
            Some(e_value) if e_value <= self.e_value_threshold => Some(alignment),
            _ => None,
        }
    }
}

#[derive(Default, Clone)]
pub struct Pipeline<S, C, A>
where
    S: SeedStep,
    C: CloudSearchStep,
    A: AlignStep,
{
    pub seed: S,
    pub cloud_search: C,
    pub align: A,
}

impl<S, C, A> Pipeline<S, C, A>
where
    S: SeedStep,
    C: CloudSearchStep,
    A: AlignStep,
{
    fn run(&mut self, profile: &mut Profile, target: &Sequence) -> Option<Alignment> {
        self.seed
            .run(profile, target)
            .and_then(|seed| self.cloud_search.run(profile, target, &seed))
            .and_then(|bounds| self.align.run(profile, target, bounds))
    }
}

pub enum HeaderStatus {
    Unwritten,
    Written,
}

pub struct Output {
    alignment_writer: Box<dyn Write + Send>,
    table_writer: Box<dyn Write + Send>,
    table_format: TableFormat,
    header_status: HeaderStatus,
}

impl Output {
    pub fn new(args: &AlignOutputArgs) -> anyhow::Result<Self> {
        Ok(Self {
            alignment_writer: Box::new(stdout()),
            table_writer: Box::new(args.tsv_results_path.open(true)?),
            table_format: TableFormat::new(&DEFAULT_COLUMNS)?,
            header_status: HeaderStatus::Unwritten,
        })
    }

    pub fn write(&mut self, alignments: &mut [Alignment]) -> anyhow::Result<()> {
        self.table_format.reset_widths();
        self.table_format.update_widths(alignments);

        alignments.sort_by(|a, b| a.e_value.partial_cmp(&b.e_value).unwrap());

        if let HeaderStatus::Unwritten = self.header_status {
            let header = TableFormat::header(&self.table_format)?;
            writeln!(self.table_writer, "{header}")?;
            self.header_status = HeaderStatus::Written;
        }

        alignments.iter().for_each(|ali| {
            writeln!(
                self.table_writer,
                "{}",
                ali.tab_string_formatted(&self.table_format)
            )
            .expect("failed to write tabular output");

            writeln!(self.alignment_writer, "{}", ali.ali_string())
                .expect("failed to write alignment output");
        });

        Ok(())
    }
}

pub fn align<S, C, A>(
    profiles: &dyn Database<Arc<Mutex<Profile>>>,
    targets: &(dyn Database<Arc<Sequence>> + Send + Sync),
    pipeline: Pipeline<S, C, A>,
    output: Arc<Mutex<Output>>,
) where
    S: SeedStep + Sync + Send,
    C: CloudSearchStep + Sync + Send,
    A: AlignStep + Sync + Send,
{
    profiles.into_par_iter().for_each_with(
        (pipeline, targets, output),
        |(pipeline, targets, output), profile| {
            let mut profile_guard = profile.lock().unwrap();
            let mut alignments: Vec<Alignment> = targets
                .iter()
                .filter_map(|target| pipeline.run(&mut profile_guard, &target))
                .collect();

            match output.lock() {
                Ok(mut guard) => {
                    guard.write(&mut alignments).unwrap();
                }
                Err(_) => panic!(),
            }
        },
    )
}
