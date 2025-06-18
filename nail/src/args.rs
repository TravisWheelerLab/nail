use std::path::PathBuf;

use clap::{Args, Parser, Subcommand};

#[derive(Subcommand)]
pub enum NailSubCommands {
    #[command(about = "Run nail's protein search pipeline")]
    Search(SearchArgs),
}

#[derive(Parser)]
#[command(version)]
#[command(name = "nail")]
#[command(
    about = "Using MMseqs2 to find rough alignment seeds, perform bounded profile HMM sequence alignment"
)]
pub struct NailCli {
    #[command(subcommand)]
    pub command: NailSubCommands,
}

#[derive(Debug, Args)]
pub struct SearchArgs {
    /// The query database file
    #[arg(value_name = "QUERY.[fasta:hmm]")]
    pub query_path: PathBuf,

    /// The target database file
    #[arg(value_name = "TARGET.fasta")]
    pub target_path: PathBuf,

    /// The number of threads that nail will use
    #[arg(short = 't', default_value_t = 8usize, value_name = "N")]
    pub num_threads: usize,

    /// Print out pipeline summary statistics
    #[arg(short = 's', action)]
    pub print_summary_stats: bool,

    /// Don't write any tabular results, write alignments to stdout
    #[arg(short = 'x', action)]
    pub ali_to_stdout: bool,

    #[command(flatten)]
    #[clap(next_help_heading = "File I/O options")]
    pub io_args: IoArgs,

    #[command(flatten)]
    #[clap(next_help_heading = "Pipeline options")]
    pub pipeline_args: PipelineArgs,

    /// Arguments that are passed to MMseqs2
    #[command(flatten)]
    #[clap(next_help_heading = "MMseqs2 options")]
    pub mmseqs_args: MmseqsArgs,

    #[command(flatten)]
    #[clap(next_help_heading = "Expert options")]
    pub expert_args: ExpertArgs,

    #[command(flatten)]
    #[clap(next_help_heading = "Dev options")]
    pub dev_args: DevArgs,
}

#[derive(Args, Debug, Clone, Default)]
pub struct IoArgs {
    /// The file where tabular output will be written
    #[arg(long = "tbl-out", default_value = "results.tbl", value_name = "PATH")]
    pub tbl_results_path: Option<PathBuf>,

    /// The file where alignment output will be written
    #[arg(long = "ali-out", default_value = None, value_name = "PATH")]
    pub ali_results_path: Option<PathBuf>,

    /// A file containing pre-computed alignment seeds
    #[arg(long = "seeds", value_name = "PATH")]
    pub seeds_input_path: Option<PathBuf>,

    /// The file where alignment seeds will be written
    #[arg(long = "seeds-out", default_value = None, value_name = "PATH")]
    pub seeds_output_path: Option<PathBuf>,

    /// The directory where intermediate files will be placed
    #[arg(long = "tmp-dir", default_value = "tmp/", value_name = "PATH")]
    pub temp_dir_path: PathBuf,

    /// Allow nail to overwrite files
    #[arg(long = "allow-overwrite", default_value_t = false)]
    pub allow_overwrite: bool,
}

#[derive(Args, Debug, Clone, Default)]
pub struct PipelineArgs {
    /// Pruning parameter alpha
    #[arg(
        short = 'A',
        default_value_t = 10.0,
        value_name = "X",
        help = "Cloud search parameter α:\n  \
                local score pruning threshold"
    )]
    pub alpha: f32,

    /// Pruning parameter beta
    #[arg(
        short = 'B',
        default_value_t = 16.0,
        value_name = "X",
        help = "Cloud search parameter β:\n  \
                global score pruning threshold"
    )]
    pub beta: f32,

    /// Pruning parameter gamma
    #[arg(
        short = 'G',
        default_value_t = 5,
        value_name = "N",
        help = "Cloud search parameter γ:\n  \
                at minimum, compute N anti-diagonals"
    )]
    pub gamma: usize,

    /// Seeding filter threshold
    #[arg(
        short = 'S',
        default_value_t = 0.01f64,
        value_name = "X",
        help = "Seeding filter threshold:\n  \
                filter hits with P-value > X"
    )]
    pub seed_pvalue_threshold: f64,

    /// Cloud search filter threshold
    #[arg(
        short = 'C',
        default_value_t = 1e-3,
        value_name = "X",
        help = "Cloud search threshold:\n  \
                filter hits with P-value > X"
    )]
    pub cloud_pvalue_threshold: f64,

    /// Forward filter threshold
    #[arg(
        short = 'F',
        default_value_t = 1e-4,
        value_name = "X",
        help = "Forward filter threshold:\n  \
                filter hits with P-value > X"
    )]
    pub forward_pvalue_threshold: f64,

    /// Final E-value threshold
    #[arg(
        short = 'E',
        default_value_t = 10.0,
        value_name = "X",
        help = "Final reporting threshold:\n  \
                filter hits with E-value > X"
    )]
    pub e_value_threshold: f64,

    /// Seed alignments twice (high/low expected sequence divergence)
    #[arg(long = "double-seed", action)]
    pub double_seed: bool,

    /// Produce alignment seeds and terminate
    #[arg(long = "only-seed", action)]
    pub only_seed: bool,
}

#[derive(Args, Debug, Clone, Default)]
pub struct ExpertArgs {
    /// Override the number of comparisons used for E-value calculation
    #[arg(short = 'Z', value_name = "N")]
    pub target_database_size: Option<usize>,

    /// Don't compute sequence composition bias score correction
    #[arg(long = "no-null2", action)]
    pub no_null_two: bool,
}

#[derive(Args, Debug, Clone, Default)]
pub struct DevArgs {
    /// Where to place stats output
    #[arg(long, value_name = "PATH", hide = true)]
    pub stats_results_path: Option<PathBuf>,

    /// Compute the full DP matrices
    #[arg(long, action, hide = true)]
    pub full_dp: bool,

    /// Give a full seed for every possible alignment
    #[arg(long, action, hide = true)]
    pub max_seed: bool,
}

#[derive(Args, Debug, Clone, Default)]
pub struct MmseqsArgs {
    /// MMseqs2 prefilter: k-mer length (0: automatically set to optimum)
    #[arg(long = "mmseqs-k", default_value_t = 0usize, value_name = "N")]
    pub k: usize,

    /// MMseqs2 prefilter: k-mer threshold for generating similar k-mer lists
    #[arg(long = "mmseqs-k-score", default_value_t = 80usize, value_name = "N")]
    pub k_score: usize,

    /// MMseqs2 prefilter: Accept only matches with ungapped alignment score above threshold
    #[arg(
        long = "mmseqs-min-ungapped-score",
        default_value_t = 15usize,
        value_name = "N"
    )]
    pub min_ungapped_score: usize,

    /// MMseqs2 prefilter: Maximum results per query sequence allowed to pass the prefilter
    #[arg(
        long = "mmseqs-max-seqs",
        default_value_t = 1000usize,
        value_name = "N"
    )]
    pub max_seqs: usize,
}
