use libnail::structs::Profile;

use crate::{
    args::SearchArgs,
    io::{Fasta, P7Hmm, SeedList, Seeds},
    mmseqs::{
        run_mmseqs_search, write_mmseqs_profile_database, write_mmseqs_sequence_database,
        MmseqsDbPaths, MmseqsScoreModel,
    },
};

pub fn seed_profile_to_sequence(
    profiles: &P7Hmm,
    seqs: &Fasta,
    args: &SearchArgs,
) -> anyhow::Result<Seeds> {
    let paths = MmseqsDbPaths::new(&args.io_args.temp_dir_path);

    write_mmseqs_sequence_database(seqs, &paths.target_db)?;

    let align_tsv = match &args.io_args.seeds_output_path {
        Some(path) => path,
        None => &paths.dir()?.join("seeds.tsv"),
    };

    write_mmseqs_profile_database(profiles, &paths.query_db)?;
    run_mmseqs_search(&paths, align_tsv, args, MmseqsScoreModel::Profile)?;
    let seeds = Seeds::from_path(align_tsv)?;
    Ok(seeds)
}

pub fn seed_sequence_to_sequence(
    queries: &Fasta,
    targets: &Fasta,
    args: &SearchArgs,
) -> anyhow::Result<Seeds> {
    let paths = MmseqsDbPaths::new(&args.io_args.temp_dir_path);

    write_mmseqs_sequence_database(targets, &paths.target_db)?;
    write_mmseqs_sequence_database(queries, &paths.query_db)?;

    let align_tsv = match &args.io_args.seeds_output_path {
        Some(path) => path,
        None => &paths.dir()?.join("seeds.tsv"),
    };

    run_mmseqs_search(&paths, align_tsv, args, MmseqsScoreModel::Blosum62)?;
    let seeds = Seeds::from_path(align_tsv)?;
    Ok(seeds)
}

dyn_clone::clone_trait_object!(SeedStage);
pub trait SeedStage: dyn_clone::DynClone + Send + Sync {
    fn run(&mut self, profile: &Profile) -> Option<SeedList>;
}

#[derive(Clone)]
pub struct DefaultSeedStage {
    seeds: Seeds,
}

impl DefaultSeedStage {
    pub fn new(seeds: Seeds) -> Self {
        DefaultSeedStage { seeds }
    }
}

impl SeedStage for DefaultSeedStage {
    fn run(&mut self, profile: &Profile) -> Option<SeedList> {
        self.seeds.get(&profile.name)
    }
}
