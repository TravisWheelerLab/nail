use crate::pipeline::{AlignArgs, SearchArgs};
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
    #[command(about = "Run the entire nail pipeline: prep, seed, & align")]
    Search(SearchArgs),
    #[command(about = "Search with the query against the target, using alignment seeds")]
    Align(AlignArgs),
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
