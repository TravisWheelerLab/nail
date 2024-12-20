use std::path::PathBuf;

use clap::{Args, Parser, Subcommand};

#[derive(Subcommand)]
pub enum SubCommands {
    #[command(about = "")]
    Search(SearchArgs),
    #[command(about = "")]
    Seed(SeedArgs),
}

#[derive(Parser)]
#[command(name = "nail")]
#[command(
    about = "Using MMseqs2 to find rough alignment seeds, perform bounded profile HMM sequence alignment"
)]
pub struct Cli {
    #[command(subcommand)]
    pub command: SubCommands,
}

#[derive(Args, Debug, Clone)]
pub struct CommonArgs {
    /// The number of threads that nail will use
    #[arg(
        short = 't',
        long = "threads",
        default_value_t = 8usize,
        value_name = "n"
    )]
    pub num_threads: usize,

    /// Allow nail to overwrite files
    #[arg(short = 'q', long = "allow-overwrite", default_value_t = false)]
    pub allow_overwrite: bool,
}

#[derive(Debug, Args)]
pub struct SearchArgs {
    /// Query file
    #[arg(value_name = "QUERY.[fasta:hmm]")]
    pub query_path: PathBuf,

    /// Target file
    #[arg(value_name = "TARGET.fasta")]
    pub target_path: PathBuf,

    /// The path to pre-computed alignment seeds
    #[arg(short = 's', long = "seeds")]
    pub seeds_path: Option<PathBuf>,

    /// Arguments that control output options
    #[command(flatten)]
    pub output_args: OutputArgs,

    /// Arguments that are passed to libnail functions
    #[command(flatten)]
    pub nail_args: NailArgs,

    /// Arguments that are passed to MMseqs2
    #[command(flatten)]
    pub mmseqs_args: MmseqsArgs,

    /// Arguments that are common across all nail subcommands
    #[command(flatten)]
    pub common_args: CommonArgs,
}

#[derive(Args, Debug, Clone)]
pub struct SeedArgs {
    /// Query file
    #[arg(value_name = "QUERY.[fasta:hmm]")]
    pub query_path: PathBuf,

    /// Target file
    #[arg(value_name = "TARGET.fasta")]
    pub target_path: PathBuf,

    /// Where to place the seeds output file
    #[arg(short = 's', long = "seeds", default_value = "seeds.json")]
    pub seeds_path: PathBuf,

    /// Arguments that are passed to MMseqs2
    #[command(flatten)]
    pub mmseqs_args: MmseqsArgs,

    /// Arguments that are common across all nail subcommands
    #[command(flatten)]
    pub common_args: CommonArgs,
}

#[derive(Args, Debug, Clone, Default)]
pub struct MmseqsArgs {
    /// The directory where intermediate files will be placed
    #[arg(long = "prep", value_name = "PATH", default_value = "prep/")]
    pub prep_dir: PathBuf,

    /// MMseqs2 prefilter: k-mer length (0: automatically set to optimum)
    #[arg(long = "mmseqs-k", default_value_t = 0usize)]
    pub k: usize,

    /// MMseqs2 prefilter: k-mer threshold for generating similar k-mer lists
    #[arg(long = "mmseqs-k-score", default_value_t = 80usize)]
    pub k_score: usize,

    /// MMseqs2 prefilter: Accept only matches with ungapped alignment score above threshold
    #[arg(long = "mmseqs-min-ungapped-score", default_value_t = 15usize)]
    pub min_ungapped_score: usize,

    /// MMseqs2 prefilter: Maximum results per query sequence allowed to pass the prefilter
    #[arg(long = "mmseqs-max-seqs", default_value_t = 1000usize)]
    pub max_seqs: usize,

    /// MMseqs2 align: Include matches below this P-value as seeds.
    ///
    // Note: the MMseqs2 align tool only allows thresholding by E-value, so the P-value supplied
    // here is multiplied by the size of the target database (i.e. number of sequences) to achieve
    // an E-value threshold that is effectively the same as the chosen P-value threshold.
    #[arg(long = "mmseqs-pvalue-threshold", default_value_t = 0.01f64)]
    pub pvalue_threshold: f64,
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
    #[arg(
        short = 'C',
        long = "cloud-thresh",
        default_value_t = 1e-3,
        value_name = "F"
    )]
    pub cloud_pvalue_threshold: f64,

    /// The P-value threshold for promoting hits past forward
    #[arg(
        short = 'F',
        long = "forward-thresh",
        default_value_t = 1e-4,
        value_name = "F"
    )]
    pub forward_pvalue_threshold: f64,

    /// Compute the full dynamic programming matrices during alignment
    #[arg(long, action)]
    pub full_dp: bool,
}

#[derive(Args, Debug, Clone)]
pub struct OutputArgs {
    /// Only report hits with an E-value below this value
    #[arg(short = 'E', default_value_t = 10.0, value_name = "F")]
    pub e_value_threshold: f64,

    /// Where to place tabular output
    #[arg(
        short = 'T',
        long = "tab-output",
        default_value = "results.tbl",
        value_name = "path"
    )]
    pub tbl_results_path: PathBuf,

    /// Where to place alignment output
    #[arg(short = 'O', long = "output", value_name = "path")]
    pub ali_results_path: Option<PathBuf>,

    /// Where to place stats output
    #[arg(short = 'S', long = "stats-output", value_name = "path")]
    pub stats_results_path: Option<PathBuf>,
}
