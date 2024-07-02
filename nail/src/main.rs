mod args;
mod mmseqs;
mod pipeline;
mod search;
mod util;

use args::{Cli, SubCommands};
use search::{search, seed};
use util::{check_mmseqs_installed, set_threads};

use clap::Parser;

fn main() -> anyhow::Result<()> {
    match Cli::parse().command {
        SubCommands::Search(args) => {
            check_mmseqs_installed()?;
            set_threads(args.common_args.num_threads)?;
            search(args)?;
        }
        SubCommands::Seed(args) => {
            check_mmseqs_installed()?;
            set_threads(args.common_args.num_threads)?;
            seed(args)?;
        }
    }
    Ok(())
}
