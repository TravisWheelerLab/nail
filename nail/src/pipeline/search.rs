use crate::args::{guess_query_format_from_query_file, FileFormat};
use crate::cli::CommonArgs;
use crate::database::{Database, ProfileCollection, SequenceCollection};
use crate::extension_traits::PathBufExt;
use crate::pipeline::{
    seed, AlignArgs, AlignOutputArgs, MmseqsArgs, NailArgs, PrepDirArgs, SeedArgs,
};
use anyhow::Context;
use clap::Args;
use libnail::structs::hmm::parse_hmms_from_p7hmm_file;
use libnail::structs::{Profile, Sequence};
use std::path::PathBuf;
use std::sync::{Arc, Mutex};

use super::{
    align, map_p7_to_mmseqs_profiles, DefaultAlignStep, DefaultCloudSearchStep, DefaultSeedStep,
    Output, Pipeline,
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

    /// Provides paths to the prep dir and the files placed in it
    #[command(flatten)]
    pub prep_dir: PrepDirArgs,

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

pub fn search(args: &SearchArgs) -> anyhow::Result<()> {
    {
        // quickly make sure we can write the results
        args.output_args.tsv_results_path.open(true)?;
        if let Some(path) = &args.output_args.ali_results_path {
            path.open(true)?;
        }
    }

    let seeds_path = args.prep_dir.path.join("./seeds.json");

    let seed_args = SeedArgs {
        prep_dir_path: Default::default(),
        seeds_path: seeds_path.clone(),
        query_path: args.query_path.clone(),
        target_path: args.target_path.clone(),
        prep_dir: args.prep_dir.clone(),
        mmseqs_args: args.mmseqs_args.clone(),
        common_args: args.common_args.clone(),
    };

    let align_args = AlignArgs {
        query_path: args.query_path.clone(),
        target_path: args.target_path.clone(),
        seeds_path,
        nail_args: args.nail_args.clone(),
        output_args: args.output_args.clone(),
        common_args: args.common_args.clone(),
    };

    let mut seeds = seed(&seed_args)?;

    let query_format = guess_query_format_from_query_file(&args.query_path)?;
    let queries: Box<dyn Database<Arc<Mutex<Profile>>>> = match query_format {
        FileFormat::Fasta => {
            let queries = SequenceCollection::new(
                Sequence::amino_from_fasta(&align_args.query_path)
                    .context("failed to read query fasta")?,
            );
            Box::new(queries)
        }
        FileFormat::Hmm => {
            let profiles = parse_hmms_from_p7hmm_file(&args.query_path)
                .context("failed to read query hmm")?
                .iter()
                .map(Profile::new)
                .collect();
            let queries = ProfileCollection::new(profiles);

            Box::new(queries)
        }
        FileFormat::Stockholm => {
            let profiles = parse_hmms_from_p7hmm_file(args.prep_dir.prep_query_hmm_path())
                .context("failed to read query hmm")?
                .iter()
                .map(Profile::new)
                .collect();
            let queries = ProfileCollection::new(profiles);

            let p7_to_mmseqs_map = map_p7_to_mmseqs_profiles(&queries, &seed_args)?;

            seeds.iter_mut().for_each(|(a, b)| {
                let map = p7_to_mmseqs_map.get(a).unwrap();
                b.values_mut().for_each(|seed| {
                    seed.profile_start = map[seed.profile_start].max(1);
                    seed.profile_end = map[seed.profile_end];
                });
            });

            Box::new(queries)
        }
        _ => {
            panic!()
        }
    };

    let targets = SequenceCollection::new(
        Sequence::amino_from_fasta(&align_args.target_path)
            .context("failed to read target fasta")?,
    );

    let pipeline = Pipeline {
        seed: DefaultSeedStep::new(seeds),
        cloud_search: DefaultCloudSearchStep::new(&align_args),
        align: DefaultAlignStep::new(&align_args, targets.len()),
    };

    let output = Arc::new(Mutex::new(Output::new(&args.output_args)?));

    align(&*queries, &targets, pipeline, output);

    Ok(())
}
