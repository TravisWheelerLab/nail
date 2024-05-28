use std::collections::HashMap;
use std::path::PathBuf;

use libnail::align::structs::Seed;

use crate::args::CommonArgs;
use crate::mmseqs::MmseqsArgs;

use clap::Args;

pub type SeedMap = HashMap<String, HashMap<String, Seed>>;

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
