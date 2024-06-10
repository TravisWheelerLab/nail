use std::path::PathBuf;

use crate::mmseqs::MmseqsArgs;
use clap::{Args, Parser, Subcommand};

#[derive(Args, Debug, Clone)]
pub struct CommonArgs {
    /// The number of threads to use
    #[arg(
        short = 't',
        long = "threads",
        default_value_t = 8usize,
        value_name = "n"
    )]
    pub num_threads: usize,
}

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

#[derive(Debug, Args)]
pub struct SearchArgs {
    /// Query file
    #[arg(value_name = "QUERY.[fasta:sto]")]
    pub query_path: PathBuf,

    /// Target file
    #[arg(value_name = "TARGET.fasta")]
    pub target_path: PathBuf,

    /// The path to a pre-built P7HMM file
    #[arg(short = 'q', long = "query-hmm", value_name = "QUERY.hmm")]
    pub prebuilt_query_hmm_path: Option<PathBuf>,

    /// Arguments that control output options
    #[command(flatten)]
    pub output_args: AlignOutputArgs,

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

    /// Arguments that are passed to MMseqs2
    #[command(flatten)]
    pub mmseqs_args: MmseqsArgs,

    /// Arguments that are common across all nail subcommands
    #[command(flatten)]
    pub common_args: CommonArgs,
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
        default_value = "results.tbl",
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
