mod args;
mod cli;
mod extension_traits;
mod pipeline;

use cli::Cli;
use extension_traits::CommandExt;
use pipeline::{align, prep, search, seed};

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

fn main() -> Result<()> {
    check_hmmer_installed()?;
    check_mmseqs_installed()?;

    match Cli::parse().command {
        SubCommands::Search(args) => {
            search(&args)?;
        }
        SubCommands::Prep(args) => {
            prep(&args)?;
        }
        SubCommands::Seed(mut args) => {
            // TODO: I'd like to think of a way to remove this nonsense
            args.prep_dir.path = args.prep_dir_path.clone();

            seed(&args)?;
        }
        SubCommands::Align(args) => {
            align(&args, None, None)?;
        }
    }
    Ok(())
}
