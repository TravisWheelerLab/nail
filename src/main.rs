use anyhow::Result;
use clap::Parser;
use nale::pipelines::{pipeline_bounded, pipeline_naive};
use nale::structs::hmm::parse_hmms_from_p7hmm_file;
use nale::structs::{Profile, Sequence};

#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    /// query hmm
    query: String,
    /// target fasta
    target: String,
}

fn main() -> Result<()> {
    let args = Args::parse();

    // println!("{:?}", args);

    let hmms = parse_hmms_from_p7hmm_file(args.query)?;
    let mut profiles: Vec<Profile> = hmms.iter().map(|hmm| Profile::new(hmm)).collect();
    let targets = Sequence::amino_from_fasta(&args.target)?;

    // pipeline_naive(&mut profiles, &targets)?;
    pipeline_bounded(&mut profiles, &targets)?;

    Ok(())
}
