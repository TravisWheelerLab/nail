use std::collections::HashMap;

use libnail::{align::structs::Seed, structs::Profile};

use crate::{
    args::SearchArgs,
    io::Fasta,
    mmseqs::{
        run_mmseqs_search, seeds_from_mmseqs_align_tsv, write_mmseqs_profile_database,
        write_mmseqs_sequence_database, MmseqsDbPaths, MmseqsScoreModel,
    },
};

fn merge_seed_maps<'a>(
    mut seed_map_a: SeedMap,
    mut seed_map_b: SeedMap,
    query_names: impl Iterator<Item = &'a str>,
) -> SeedMap {
    query_names.for_each(|query: &str| {
        let seeds_a = seed_map_a.get_mut(query);
        let seeds_b = seed_map_b.remove(query);
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
                seed_map_a.insert(query.to_string(), b);
            }
            _ => {}
        }
    });

    seed_map_a
}

pub fn seed_profile_to_sequence(
    queries: &[Profile],
    targets: &Fasta,
    args: &SearchArgs,
) -> anyhow::Result<SeedMap> {
    let paths = MmseqsDbPaths::new(&args.io_args.temp_dir_path);

    write_mmseqs_sequence_database(targets, &paths.target_db)?;
    write_mmseqs_profile_database(queries, &paths.query_db)?;

    run_mmseqs_search(&paths, args, MmseqsScoreModel::Profile)?;

    let seed_map_a = seeds_from_mmseqs_align_tsv(&paths.align_tsv)?;

    if !args.pipeline_args.double_seed {
        return Ok(seed_map_a);
    }

    let queries_b: Vec<_> = queries
        .iter()
        .filter(|p| p.relative_entropy() < 1.0)
        .map(|p| {
            let mut p2 = p.clone();
            p2.adjust_mean_relative_entropy(1.0).unwrap();
            p2
        })
        .collect();

    write_mmseqs_profile_database(&queries_b, &paths.query_db)?;

    run_mmseqs_search(&paths, args, MmseqsScoreModel::Profile)?;

    let seed_map_b = seeds_from_mmseqs_align_tsv(&paths.align_tsv)?;

    let names = queries.iter().map(|q| q.name.as_str());
    let seed_map_merged = merge_seed_maps(seed_map_a, seed_map_b, names);

    Ok(seed_map_merged)
}

pub fn seed_sequence_to_sequence(
    queries: &Fasta,
    targets: &Fasta,
    args: &SearchArgs,
) -> anyhow::Result<SeedMap> {
    let paths = MmseqsDbPaths::new(&args.io_args.temp_dir_path);

    write_mmseqs_sequence_database(targets, &paths.target_db)?;
    write_mmseqs_sequence_database(queries, &paths.query_db)?;

    run_mmseqs_search(&paths, args, MmseqsScoreModel::Blosum62)?;

    let seed_map_a = seeds_from_mmseqs_align_tsv(&paths.align_tsv)?;

    if !args.pipeline_args.double_seed {
        return Ok(seed_map_a);
    }

    run_mmseqs_search(&paths, args, MmseqsScoreModel::Blosum80)?;

    let seed_map_b = seeds_from_mmseqs_align_tsv(&paths.align_tsv)?;

    let names = queries.names_iter();
    let seed_map_merged = merge_seed_maps(seed_map_a, seed_map_b, names);

    Ok(seed_map_merged)
}

dyn_clone::clone_trait_object!(SeedStage);
pub trait SeedStage: dyn_clone::DynClone + Send + Sync {
    fn run(&mut self, profile: &Profile) -> Option<&HashMap<String, Seed>>;
}

pub type SeedMap = HashMap<String, HashMap<String, Seed>>;

#[derive(Default, Clone)]
pub struct DefaultSeedStage {
    seeds: SeedMap,
}

impl DefaultSeedStage {
    pub fn new(seeds: SeedMap) -> Self {
        DefaultSeedStage { seeds }
    }
}

impl SeedStage for DefaultSeedStage {
    fn run(&mut self, profile: &Profile) -> Option<&HashMap<String, Seed>> {
        self.seeds.get(&profile.name)
    }
}
