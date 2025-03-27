mod args;
mod io;
mod mmseqs;
mod pipeline;
mod search;
mod stats;
mod util;

use args::{NailCli, NailSubCommands};
use search::search;
use util::{check_mmseqs_installed, set_threads};

use clap::Parser;

#[cfg(feature = "jemalloc")]
use jemallocator::Jemalloc;

#[cfg(feature = "jemalloc")]
#[global_allocator]
static GLOBAL: Jemalloc = Jemalloc;

fn main() -> anyhow::Result<()> {
    color_backtrace::install();
    match NailCli::parse().command {
        NailSubCommands::Search(args) => {
            check_mmseqs_installed()?;
            set_threads(args.num_threads)?;
            search(args)?;
        }
    }

    Ok(())
}
