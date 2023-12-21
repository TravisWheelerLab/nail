use crate::cli::CommonArgs;
use crate::extension_traits::PathBufExt;
use crate::pipeline::{
    align, prep, seed, AlignArgs, AlignOutputArgs, MmseqsArgs, NailArgs, PrepArgs, PrepDirArgs,
    SeedArgs,
};
use clap::Args;
use std::path::PathBuf;

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
    let (profiles, seed_map) = seed(&seed_args)?;

    align(&align_args, Some(profiles), Some(seed_map))?;
    Ok(())
}
