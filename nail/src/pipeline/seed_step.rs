use std::{collections::HashMap, path::Path};

use libnail::{
    align::structs::Seed,
    structs::{Profile, Sequence},
};

use crate::{
    mmseqs::{consts::ALIGN_TSV, run_mmseqs_search, seeds_from_mmseqs_align_tsv, MmseqsArgs},
    search::Queries,
};

pub trait SeedStep: dyn_clone::DynClone {
    fn run(
        &self,
        profile: &Profile,
        target: &HashMap<String, Sequence>,
    ) -> Option<&HashMap<String, Seed>>;
}

dyn_clone::clone_trait_object!(SeedStep);

pub type SeedMap = HashMap<String, HashMap<String, Seed>>;

#[derive(Default, Clone)]
pub struct DefaultSeedStep {
    seeds: SeedMap,
}

impl DefaultSeedStep {
    pub fn new(
        queries: &Queries,
        targets: &[Sequence],
        prep_dir: &impl AsRef<Path>,
        num_threads: usize,
        mmseqs_args: &MmseqsArgs,
    ) -> anyhow::Result<Self> {
        run_mmseqs_search(queries, targets, prep_dir, num_threads, mmseqs_args)?;

        let seeds = seeds_from_mmseqs_align_tsv(prep_dir.as_ref().join(ALIGN_TSV))?;

        Ok(Self { seeds })
    }

    pub fn from_seed_map(seeds: SeedMap) -> Self {
        DefaultSeedStep { seeds }
    }
}

impl SeedStep for DefaultSeedStep {
    fn run(
        &self,
        profile: &Profile,
        _target: &HashMap<String, Sequence>,
    ) -> Option<&HashMap<String, Seed>> {
        self.seeds.get(&profile.name)
    }
}
