use anyhow::Result;
use clap::ArgAction;
use clap::Parser;
use nale::align::bounded::structs::{CloudBoundGroup, CloudSearchParams};
use nale::output::output_tabular::write_tabular_output;
use nale::pipelines::{
    pipeline_bounded, pipeline_naive, BoundedPipelineParams, NaivePipelineParams,
};
use nale::structs::hmm::parse_hmms_from_p7hmm_file;
use nale::structs::{Profile, Sequence};
use nale::viz::SodaJson;
use std::io;
use std::io::stdout;
use std::path::PathBuf;

#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
pub struct Args {
    /// Query hmm file
    query: String,
    /// Target fasta file
    target: String,
    /// Path for alignment output
    #[arg(long)]
    ali_out: Option<String>,
    /// Path for tabular output
    table_out: Option<String>,
    /// Write debugging information
    #[arg(long, action = ArgAction::SetTrue)]
    debug: Option<bool>,
    /// Allow output files to be overwritten
    #[arg(long, action = ArgAction::SetTrue)]
    allow_overwrite: Option<bool>,
}

impl Args {
    pub fn params_naive(&self) -> NaivePipelineParams {
        NaivePipelineParams {
            write_debug: match self.debug {
                Some(val) => val,
                None => false,
            },
            allow_overwrite: match self.allow_overwrite {
                Some(val) => val,
                None => false,
            },
            // TODO: parameterize this
            root_debug_dir_path: PathBuf::from("./nale-debug"),
        }
    }

    pub fn params_bounded(&self) -> BoundedPipelineParams {
        BoundedPipelineParams {
            write_debug: match self.debug {
                Some(val) => val,
                None => false,
            },
            allow_overwrite: match self.allow_overwrite {
                Some(val) => val,
                None => false,
            },
            // TODO: parameterize this
            root_debug_dir_path: PathBuf::from("./nale-debug"),
            cloud_search_params: CloudSearchParams {
                target_start: 0,
                target_end: 0,
                profile_start: 0,
                profile_end: 0,
                gamma: 0,
                alpha: 0.0,
                beta: 0.0,
            },
        }
    }
}

fn main() -> Result<()> {
    let mut b1 = CloudBoundGroup::new(100, 100);
    let mut b2 = CloudBoundGroup::new(100, 100);

    b1.set(0, 0, 0, 0, 0);
    b1.set(1, 1, 0, 0, 1);
    b1.set(2, 2, 0, 0, 2);
    b1.set(3, 3, 0, 0, 3);
    b1.set(4, 4, 0, 0, 4);

    b2.set(200, 100, 100, 100, 100);
    b2.set(199, 100, 99, 99, 100);
    b2.set(198, 100, 98, 98, 100);
    b2.set(197, 100, 97, 97, 100);
    b2.set(196, 100, 96, 96, 100);

    b1.soda_json(&mut stdout())?;
    println!();
    b2.soda_json(&mut stdout())?;
    println!();
    CloudBoundGroup::join_bounds(&mut b1, &b2)?;
    println!();
    b1.soda_json(&mut stdout())?;
    // let args = Args::parse();
    // let hmms = parse_hmms_from_p7hmm_file(&args.query)?;
    // let mut profiles: Vec<Profile> = hmms.iter().map(|hmm| Profile::new(hmm)).collect();
    // let targets = Sequence::amino_from_fasta(&args.target)?;
    //
    // let params_naive = args.params_naive();
    // let alignments_naive = pipeline_naive(&mut profiles, &targets, &params_naive)?;
    //
    // let params_bounded = args.params_bounded();
    // let alignments_bounded = pipeline_bounded(&mut profiles, &targets, &params_bounded)?;
    //
    // write_tabular_output(&alignments_bounded, &mut io::stdout())?;
    Ok(())
}
