use anyhow::Result;
use clap::ArgAction;
use clap::Parser;
use nale::align::bounded::structs::CloudSearchParams;
use nale::output::output_tabular::write_tabular_output;
use nale::pipelines::{
    pipeline_bounded, pipeline_naive, BoundedPipelineParams, NaivePipelineParams,
};
use nale::structs::hmm::parse_hmms_from_p7hmm_file;
use nale::structs::{Profile, Sequence};
use std::io;
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
    let args = Args::parse();
    let hmms = parse_hmms_from_p7hmm_file(&args.query)?;
    let mut profiles: Vec<Profile> = hmms.iter().map(|hmm| Profile::new(hmm)).collect();
    let targets = Sequence::amino_from_fasta(&args.target)?;

    let params_naive = args.params_naive();
    let alignments_naive = pipeline_naive(&mut profiles, &targets, &params_naive)?;

    let params_bounded = args.params_bounded();
    let alignments_bounded = pipeline_bounded(&mut profiles, &targets, &params_bounded)?;

    write_tabular_output(&alignments_bounded, &mut io::stdout())?;
    Ok(())
}
