use std::{collections::HashMap, path::Path};

use anyhow::bail;
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

#[derive(Default, Clone)]
pub struct DoubleSeedStep {
    seeds: SeedMap,
}

impl DoubleSeedStep {
    pub fn new(
        queries: &Queries,
        targets: &[Sequence],
        prep_dir: &impl AsRef<Path>,
        num_threads: usize,
        mmseqs_args: &MmseqsArgs,
    ) -> anyhow::Result<Self> {
        let (queries_a, queries_b) = match queries {
            Queries::DoubleProfile(queries) => queries,
            _ => bail!("wrong query type"),
        };

        run_mmseqs_search(
            &Queries::Profile(queries_a.to_vec()),
            targets,
            prep_dir,
            num_threads,
            mmseqs_args,
        )?;
        let mut seed_map_a = seeds_from_mmseqs_align_tsv(prep_dir.as_ref().join(ALIGN_TSV))?;

        run_mmseqs_search(
            &Queries::Profile(queries_b.to_vec()),
            targets,
            prep_dir,
            num_threads,
            mmseqs_args,
        )?;
        let mut seed_map_b = seeds_from_mmseqs_align_tsv(prep_dir.as_ref().join(ALIGN_TSV))?;

        queries_a.iter().for_each(|profile| {
            let seeds_a = seed_map_a.get_mut(&profile.name);
            let seeds_b = seed_map_b.remove(&profile.name);

            match (seeds_a, seeds_b) {
                (Some(seeds_a), Some(seeds_b)) => {
                    seeds_b.into_iter().for_each(|(target, seed_b)| {
                        match seeds_a.get(&target) {
                            Some(seed_a) => seeds_a.insert(
                                target,
                                Seed {
                                    target_start: seed_a.target_start.min(seed_b.target_start),
                                    target_end: seed_a.target_end.max(seed_b.target_end),
                                    profile_start: seed_a.profile_start.min(seed_b.profile_start),
                                    profile_end: seed_a.profile_end.max(seed_b.profile_end),
                                },
                            ),
                            None => seeds_a.insert(target, seed_b),
                        };
                    });
                }
                (None, Some(b)) => {
                    seed_map_a.insert(profile.name.clone(), b);
                }
                _ => {}
            }
        });

        Ok(Self { seeds: seed_map_a })
    }
}

impl SeedStep for DoubleSeedStep {
    fn run(
        &self,
        profile: &Profile,
        _target: &HashMap<String, Sequence>,
    ) -> Option<&HashMap<String, Seed>> {
        self.seeds.get(&profile.name)
    }
}
