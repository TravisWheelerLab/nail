mod args;
mod io;
mod mmseqs;
mod pipeline;
mod search;
mod stats;
mod util;

use args::{NailCli, NailSubCommands};
use search::search;
use util::{check_mmseqs_installed, set_threads, term::*};

use clap::Parser;

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
        NailSubCommands::Search(args) => {
            args.validate()?;
            check_mmseqs_installed()?;
            set_threads(args.num_threads)?;
            search(args)?;
        }
        NailSubCommands::Dev => {}
    }

    Ok(())
}
