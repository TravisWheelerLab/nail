mod args;
mod mmseqs;
mod pipeline;
mod util;
mod viz;

use args::{Cli, SubCommands};
use pipeline::search;
use util::{check_hmmer_installed, check_mmseqs_installed, set_threads};

use clap::Parser;

fn main() -> anyhow::Result<()> {
    match Cli::parse().command {
        SubCommands::Search(args) => {
            check_hmmer_installed()?;
            check_mmseqs_installed()?;
            set_threads(args.common_args.num_threads)?;
            search(&args)?;
        }
        SubCommands::Seed(mut args) => {
            check_hmmer_installed()?;
            check_mmseqs_installed()?;

            // seed(&args)?;
        }
        SubCommands::Align(args) => {
            set_threads(args.common_args.num_threads)?;

            // align(&args, None, None)?;
        }
    }
    Ok(())
}
