mod args;
mod cli;
mod database;
mod extension_traits;
mod id;
mod pipeline;
mod viz;

use cli::Cli;
use extension_traits::CommandExt;
use pipeline::{align, search, seed};

use crate::cli::SubCommands;
use anyhow::{Context, Result};
use clap::Parser;

fn check_hmmer_installed() -> Result<()> {
    std::process::Command::new("hmmbuild")
        .arg("-h")
        .run()
        .context("hmmbuild does not appear to be in the system path")
}

fn check_mmseqs_installed() -> Result<()> {
    std::process::Command::new("mmseqs")
        .arg("-h")
        .run()
        .context("mmseqs2 does not appear to be in the system path")
}

fn set_threads(num_threads: usize) -> Result<()> {
    rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build_global()
        .context("failed to build rayon global threadpool")
}

fn main() -> Result<()> {
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
            // TODO: I'd like to think of a way to remove this nonsense
            args.prep_dir.path = args.prep_dir_path.clone();

            // seed(&args)?;
        }
        SubCommands::Align(args) => {
            set_threads(args.common_args.num_threads)?;

            // align(&args, None, None)?;
        }
    }
    Ok(())
}
