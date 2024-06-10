use std::collections::HashMap;

use libnail::{
    align::structs::Seed,
    structs::{Profile, Sequence},
};

use crate::mmseqs::{
    run_mmseqs_search, seeds_from_mmseqs_align_tsv, write_mmseqs_profile_database,
    write_mmseqs_sequence_database, MmseqsArgs,
};

fn merge_seed_maps(
    mut seed_map_a: SeedMap,
    mut seed_map_b: SeedMap,
    queries: &[impl AsRef<Profile>],
) -> SeedMap {
    queries.iter().map(|q| q.as_ref()).for_each(|profile| {
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

    seed_map_a
}

pub trait SeedStep: dyn_clone::DynClone {
    fn run(
        &mut self,
        profile: &Profile,
        target: &HashMap<String, Sequence>,
    ) -> Option<&HashMap<String, Seed>>;
}

dyn_clone::clone_trait_object!(SeedStep);

pub type SeedMap = HashMap<String, HashMap<String, Seed>>;

#[derive(Default, Clone)]
pub struct ProfileSeedStep {
    seeds: SeedMap,
}

impl ProfileSeedStep {
    pub fn new(
        queries: &[Profile],
        targets: &[Sequence],
        num_threads: usize,
        mmseqs_args: &MmseqsArgs,
    ) -> anyhow::Result<Self> {
        let query_db = mmseqs_args.prep_dir.join("queryDB");
        let target_db = mmseqs_args.prep_dir.join("targetDB");
        let align_db = mmseqs_args.prep_dir.join("alignDB");
        let align_tsv = mmseqs_args.prep_dir.join("align.tsv");

        write_mmseqs_sequence_database(targets, &target_db)?;

        write_mmseqs_profile_database(queries, &query_db)?;

        run_mmseqs_search(
            &query_db,
            &target_db,
            &align_db,
            &align_tsv,
            targets.len(),
            num_threads,
            mmseqs_args,
        )?;

        let seed_map_a = seeds_from_mmseqs_align_tsv(&align_tsv)?;

        let other_queries: Vec<_> = queries
            .iter()
            .filter(|p| p.relative_entropy() < 1.0)
            .map(|p| {
                let mut p2 = p.clone();
                p2.tune_relative_entropy(1.0);
                p2
            })
            .collect();

        write_mmseqs_profile_database(&other_queries, &query_db)?;

        run_mmseqs_search(
            query_db,
            target_db,
            align_db,
            &align_tsv,
            targets.len(),
            num_threads,
            mmseqs_args,
        )?;

        let seed_map_b = seeds_from_mmseqs_align_tsv(align_tsv)?;

        let seeds = merge_seed_maps(seed_map_a, seed_map_b, queries);

        Ok(Self { seeds })
    }

    pub fn from_seed_map(seeds: SeedMap) -> Self {
        ProfileSeedStep { seeds }
    }
}

impl SeedStep for ProfileSeedStep {
    fn run(
        &mut self,
        profile: &Profile,
        _target: &HashMap<String, Sequence>,
    ) -> Option<&HashMap<String, Seed>> {
        self.seeds.get(&profile.name)
    }
}

#[derive(Default, Clone)]
pub struct SequenceSeedStep {
    seeds: SeedMap,
}

impl SequenceSeedStep {
    pub fn new(
        queries: &[Sequence],
        targets: &[Sequence],
        num_threads: usize,
        mmseqs_args: &MmseqsArgs,
    ) -> anyhow::Result<Self> {
        let query_db = mmseqs_args.prep_dir.join("queryDB");
        let target_db = mmseqs_args.prep_dir.join("targetDB");
        let align_db = mmseqs_args.prep_dir.join("alignDB");
        let align_tsv = mmseqs_args.prep_dir.join("align.tsv");

        write_mmseqs_sequence_database(targets, &target_db)?;
        write_mmseqs_sequence_database(queries, &query_db)?;

        run_mmseqs_search(
            query_db,
            target_db,
            align_db,
            &align_tsv,
            targets.len(),
            num_threads,
            mmseqs_args,
        )?;

        let seeds = seeds_from_mmseqs_align_tsv(&align_tsv)?;

        Ok(Self { seeds })
    }

    pub fn from_seed_map(seeds: SeedMap) -> Self {
        SequenceSeedStep { seeds }
    }
}

impl SeedStep for SequenceSeedStep {
    fn run(
        &mut self,
        profile: &Profile,
        _target: &HashMap<String, Sequence>,
    ) -> Option<&HashMap<String, Seed>> {
        self.seeds.get(&profile.name)
    }
}
