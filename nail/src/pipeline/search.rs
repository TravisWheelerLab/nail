use crate::cli::CommonArgs;
use crate::database::{ProfileCollection, SequenceCollection};
use crate::extension_traits::PathBufExt;
use crate::pipeline::{
    prep, seed, AlignArgs, AlignOutputArgs, MmseqsArgs, NailArgs, PrepArgs, PrepDirArgs, SeedArgs,
};
use clap::Args;
use libnail::align::structs::Seed;
use libnail::structs::Sequence;
use std::collections::HashMap;
use std::path::PathBuf;
use std::sync::{Arc, Mutex};

use super::{align, DefaultAlignStep, DefaultCloudSearchStep, DefaultSeedStep, Output, Pipeline};

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
    }

    let seeds_path = args.prep_dir.path.join("./seeds.json");
    let prep_args = PrepArgs {
        query_path: args.query_path.clone(),
        target_path: args.target_path.clone(),
        skip_hmmbuild: args.prebuilt_query_hmm_path.is_some(),
        prep_dir: args.prep_dir.clone(),
        common_args: args.common_args.clone(),
    };

    let seed_args = SeedArgs {
        prep_dir_path: Default::default(),
        seeds_path: seeds_path.clone(),
        prebuilt_query_hmm_path: args.prebuilt_query_hmm_path.clone(),
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

    prep(&prep_args)?;
    let (profiles_vec, seed_map) = seed(&seed_args)?;

    // TODO: change seed map structure in seed()
    let mut seeds: HashMap<String, HashMap<String, Seed>> = HashMap::new();
    seed_map.iter().for_each(|(profile, seed_vec)| {
        let profile_map = seeds.entry(profile.to_string()).or_default();
        seed_vec.iter().for_each(|seed| {
            profile_map.insert(seed.target_name.clone(), seed.clone());
        });
    });

    let profiles = ProfileCollection::new(profiles_vec);
    let targets = SequenceCollection::new(Sequence::amino_from_fasta(&align_args.target_path)?);

    let pipeline = Pipeline {
        seed: DefaultSeedStep::new(seeds),
        cloud_search: DefaultCloudSearchStep::new(&align_args),
        align: DefaultAlignStep::new(&align_args, targets.len()),
    };

    let output = Arc::new(Mutex::new(Output::new(&args.output_args)?));

    align(&profiles, &targets, pipeline, output);

    Ok(())
}
