use crate::args::CommonArgs;
use crate::mmseqs::MmseqsArgs;
use crate::pipeline::{AlignArgs, AlignOutputArgs, NailArgs};
use crate::util::{guess_query_format_from_query_file, FileFormat, PathBufExt};
use anyhow::Context;
use clap::Args;
use libnail::align::structs::Seed;
use libnail::structs::hmm::parse_hmms_from_p7hmm_file;
use libnail::structs::{Profile, Sequence};
use std::collections::HashMap;
use std::path::PathBuf;
use std::sync::{Arc, Mutex};

use super::{
    run_pipeline_profile_to_sequence, run_pipeline_sequence_to_sequence, DefaultAlignStep,
    DefaultCloudSearchStep, DefaultSeedStep, FullDpCloudSearchStep, Output, Pipeline,
};

#[derive(Debug, Args)]
pub struct SearchArgs {
    /// Query file
    #[arg(value_name = "QUERY.[fasta:sto]")]
    pub query_path: PathBuf,

    /// Target file
    #[arg(value_name = "TARGET.fasta")]
    pub target_path: PathBuf,

    /// The path to a pre-built P7HMM file
    #[arg(short = 'q', long = "query-hmm", value_name = "QUERY.hmm")]
    pub prebuilt_query_hmm_path: Option<PathBuf>,

    /// The directory where intermediate files will be placed
    #[arg(long = "prep", value_name = "<PATH>", default_value = "prep/")]
    pub prep_dir: PathBuf,

    /// Arguments that control output options
    #[command(flatten)]
    pub output_args: AlignOutputArgs,

    /// Arguments that are passed to libnail functions
    #[command(flatten)]
    pub nail_args: NailArgs,

    /// Arguments that are passed to MMseqs2
    #[command(flatten)]
    pub mmseqs_args: MmseqsArgs,

    /// Arguments that are common across all nail subcommands
    #[command(flatten)]
    pub common_args: CommonArgs,
}

// TODO: move this
pub type SeedMap = HashMap<String, HashMap<String, Seed>>;

pub enum Queries {
    Sequence(Vec<Sequence>),
    Profile(Vec<Profile>),
}

pub fn search(args: &SearchArgs) -> anyhow::Result<()> {
    {
        // quickly make sure we can write the results
        args.output_args.tsv_results_path.open(true)?;
        if let Some(path) = &args.output_args.ali_results_path {
            path.open(true)?;
        }
    }

    let seeds_path = args.prep_dir.join("./seeds.json");

    let align_args = AlignArgs {
        query_path: args.query_path.clone(),
        target_path: args.target_path.clone(),
        seeds_path,
        nail_args: args.nail_args.clone(),
        output_args: args.output_args.clone(),
        common_args: args.common_args.clone(),
    };

    let query_format = guess_query_format_from_query_file(&args.query_path)?;

    let queries = match query_format {
        FileFormat::Fasta => {
            let queries = Sequence::amino_from_fasta(&align_args.query_path)
                .context("failed to read query fasta")?;

            Queries::Sequence(queries)
        }
        FileFormat::Hmm => {
            let queries: Vec<Profile> = parse_hmms_from_p7hmm_file(&args.query_path)
                .context("failed to read query hmm")?
                .iter()
                .map(Profile::new)
                .collect();
            Queries::Profile(queries)
        }
        _ => {
            panic!()
        }
    };

    let targets = Sequence::amino_from_fasta(&align_args.target_path)
        .context("failed to read target fasta")?;

    let pipeline = Pipeline {
        seed: Box::new(DefaultSeedStep::new(
            &queries,
            &targets,
            &args.prep_dir,
            args.common_args.num_threads,
            &args.mmseqs_args,
        )?),
        cloud_search: match args.nail_args.full_dp {
            true => Box::<FullDpCloudSearchStep>::default(),
            false => Box::new(DefaultCloudSearchStep::new(&align_args)),
        },
        align: Box::new(DefaultAlignStep::new(&align_args, targets.len())),
    };

    let output = Arc::new(Mutex::new(Output::new(&args.output_args)?));

    let target_map: HashMap<String, Sequence> =
        targets.into_iter().map(|t| (t.name.clone(), t)).collect();

    match queries {
        Queries::Sequence(queries) => {
            run_pipeline_sequence_to_sequence(&queries, &target_map, pipeline, output);
        }
        Queries::Profile(mut queries) => {
            run_pipeline_profile_to_sequence(&mut queries, &target_map, pipeline, output);
        }
    }

    Ok(())
}
