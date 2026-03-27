mod args;
mod dev;
mod io;
mod mmseqs;
mod pipeline;
mod search;
mod stats;
mod util;

use clap::Parser;

use crate::{
    args::{DevSubCommands, NailCli, NailSubCommands},
    dev::{dev_mx, dev_play, dev_search},
    search::search,
    util::{check_mmseqs_installed, set_threads, term::*},
};

#[cfg(feature = "jemalloc")]
#[global_allocator]
static GLOBAL: jemallocator::Jemalloc = jemallocator::Jemalloc;

fn main() {
    if let Err(e) = run() {
        eprintln!("\n{RED}error:{RESET} {e:?}\n");
        std::process::exit(1);
    }
}

fn run() -> anyhow::Result<()> {
    color_backtrace::install();
    match NailCli::parse().command {
        NailSubCommands::Search(mut args) => {
            args.validate()?;
            check_mmseqs_installed(&args.mmseqs_path)?;
            set_threads(args.num_threads)?;
            search(args)?;
        }
        NailSubCommands::Dev(cmd) => match cmd {
            DevSubCommands::Play(mut args) => {
                args.validate()?;
                check_mmseqs_installed(&args.mmseqs_path)?;
                set_threads(args.num_threads)?;
                dev_play(args)?;
            }
            DevSubCommands::Search(mut args) => {
                args.validate()?;
                check_mmseqs_installed(&args.mmseqs_path)?;
                set_threads(args.num_threads)?;
                dev_search(args)?;
            }
            DevSubCommands::Mx(mut args) => {
                args.validate()?;
                check_mmseqs_installed(&args.mmseqs_path)?;
                set_threads(args.num_threads)?;
                dev_mx(args)?;
            }
        },
    }

    Ok(())
}
