mod args;
mod io;
mod mmseqs;
mod pipeline;
mod search;
mod stats;
mod util;

use args::{NailCli, NailSubCommands};
use io::Fasta;
use search::search;
use util::{check_mmseqs_installed, set_threads};

use clap::Parser;

#[cfg(feature = "jemalloc")]
#[global_allocator]
static GLOBAL: jemallocator::Jemalloc = jemallocator::Jemalloc;

fn main() -> anyhow::Result<()> {
    color_backtrace::install();
    match NailCli::parse().command {
        NailSubCommands::Search(args) => {
            check_mmseqs_installed()?;
            set_threads(args.num_threads)?;
            search(args)?;
        }
        NailSubCommands::Dev => {
            const N: usize = 8;
            const F: &str = "./target.fa";
            const H: &str = "./query.hmm";

            let now = std::time::Instant::now();
            let a = Fasta::from_path(F)?;
            println!("{:?}", now.elapsed());
        }
    }

    Ok(())
}
