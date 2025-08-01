mod args;
mod io;
mod mmseqs;
mod pipeline;
mod search;
mod stats;
mod util;

use std::{fs::File, io::BufWriter};

use args::{NailCli, NailSubCommands};
use libnail::structs::{Hmm, Profile};
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
            let mut out = BufWriter::new(File::create(format!(
                "{}.npb",
                args.query_path.file_stem().unwrap().to_str().unwrap()
            ))?);
            Hmm::from_p7hmm(File::open(args.query_path)?)?
                .iter()
                .map(Profile::new)
                .for_each(|p| p.serialize(&mut out).unwrap());
            //
            // check_mmseqs_installed()?;
            // set_threads(args.num_threads)?;
            // search(args)?;
        }
    }

    Ok(())
}
