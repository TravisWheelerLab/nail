use std::collections::HashMap;

use libnail::{
    align::structs::Seed,
    structs::{Profile, Sequence},
};

use crate::{
    args::MmseqsArgs,
    mmseqs::{
        run_mmseqs_search, seeds_from_mmseqs_align_tsv, write_mmseqs_profile_database,
        write_mmseqs_sequence_database, MmseqsDbPaths,
    },
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
                        Some(seed_a) if seed_a.score > seed_b.score => {}
                        _ => {
                            seeds_a.insert(target, seed_b);
                        }
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

pub fn seed_profile_to_sequence(
    queries: &[Profile],
    targets: &[Sequence],
    num_threads: usize,
    mmseqs_args: &MmseqsArgs,
) -> anyhow::Result<SeedMap> {
    let paths = MmseqsDbPaths::new(&mmseqs_args.prep_dir);

    write_mmseqs_sequence_database(targets, &paths.target_db)?;
    write_mmseqs_profile_database(queries, &paths.query_db)?;

    run_mmseqs_search(&paths, targets.len(), num_threads, mmseqs_args)?;

    let seed_map_a = seeds_from_mmseqs_align_tsv(&paths.align_tsv)?;

    let queries_b: Vec<_> = queries
        .iter()
        .filter(|p| p.relative_entropy() < 1.0)
        .map(|p| {
            let mut p2 = p.clone();
            p2.raise_relative_entropy(1.0);
            p2
        })
        .collect();

    write_mmseqs_profile_database(&queries_b, &paths.query_db)?;

    run_mmseqs_search(&paths, targets.len(), num_threads, mmseqs_args)?;

    let seed_map_b = seeds_from_mmseqs_align_tsv(&paths.align_tsv)?;

    let seeds = merge_seed_maps(seed_map_a, seed_map_b, queries);

    Ok(seeds)
}

pub fn seed_sequence_to_sequence(
    queries: &[Sequence],
    targets: &[Sequence],
    num_threads: usize,
    mmseqs_args: &MmseqsArgs,
) -> anyhow::Result<SeedMap> {
    let paths = MmseqsDbPaths::new(&mmseqs_args.prep_dir);

    write_mmseqs_sequence_database(targets, &paths.target_db)?;
    write_mmseqs_sequence_database(queries, &paths.query_db)?;

    run_mmseqs_search(&paths, targets.len(), num_threads, mmseqs_args)?;

    let seeds = seeds_from_mmseqs_align_tsv(&paths.align_tsv)?;

    Ok(seeds)
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
pub struct DefaultSeedStep {
    seeds: SeedMap,
}

impl DefaultSeedStep {
    pub fn new(seeds: SeedMap) -> Self {
        DefaultSeedStep { seeds }
    }
}

impl SeedStep for DefaultSeedStep {
    fn run(
        &mut self,
        profile: &Profile,
        _target: &HashMap<String, Sequence>,
    ) -> Option<&HashMap<String, Seed>> {
        self.seeds.get(&profile.name)
    }
}
