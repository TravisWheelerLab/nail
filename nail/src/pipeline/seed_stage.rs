use std::{
    collections::HashMap,
    io::{BufRead, BufReader, BufWriter, Read, Write},
};

use libnail::{align::structs::Seed, structs::Profile};

use crate::{
    args::SearchArgs,
    io::{Fasta, P7Hmm, ProfileDatabase, SequenceDatabase},
    mmseqs::{
        run_mmseqs_search, seeds_from_mmseqs_align_tsv, write_mmseqs_profile_database,
        write_mmseqs_sequence_database, MmseqsDbPaths, MmseqsScoreModel,
    },
    search::Queries,
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
    profiles: &P7Hmm,
    seqs: &Fasta,
    args: &SearchArgs,
) -> anyhow::Result<SeedMap> {
    let paths = MmseqsDbPaths::new(&args.io_args.temp_dir_path);

    write_mmseqs_sequence_database(seqs, &paths.target_db)?;
    write_mmseqs_profile_database(profiles, &paths.query_db, None)?;

    run_mmseqs_search(&paths, args, MmseqsScoreModel::Profile)?;
    let seed_map_a = seeds_from_mmseqs_align_tsv(&paths.align_tsv)?;

    if !args.pipeline_args.double_seed {
        return Ok(seed_map_a);
    }

    write_mmseqs_profile_database(&profiles, &paths.query_db, Some(1.0))?;

    run_mmseqs_search(&paths, args, MmseqsScoreModel::Profile)?;
    let seed_map_b = seeds_from_mmseqs_align_tsv(&paths.align_tsv)?;

    let seed_map_merged = merge_seed_maps(seed_map_a, seed_map_b, profiles.names_iter());

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

pub fn write_seed_map(map: &SeedMap, buf: &mut impl Write) -> anyhow::Result<()> {
    let mut writer = BufWriter::new(buf);
    map.iter().try_for_each(|(prf_name, seeds)| {
        writeln!(writer, "{}", prf_name)?;
        seeds.iter().try_for_each(|(seq_name, seed)| {
            writeln!(writer, ">{}", seq_name)?;
            writeln!(writer, "{}", seed)
        })?;
        writeln!(writer, "//")
    })?;

    Ok(())
}

pub fn read_seed_map(buf: &mut impl Read) -> anyhow::Result<SeedMap> {
    let mut seeds: SeedMap = HashMap::new();
    let mut reader = BufReader::new(buf);
    let mut line = String::new();

    enum ParseState {
        PrfName,
        SeqName,
        Seed,
        End,
    }

    let mut state = ParseState::PrfName;
    let mut prf_name = "".to_string();
    let mut seq_name = "".to_string();
    let mut prf_seeds = HashMap::new();
    while reader.read_line(&mut line).unwrap() > 0 {
        let trimmed = line.trim_end();
        if trimmed.is_empty() {
            continue;
        } else if trimmed == "//" {
            state = ParseState::End;
        }

        match state {
            ParseState::PrfName => {
                prf_name = trimmed.to_string();
                state = ParseState::SeqName;
            }
            ParseState::SeqName => {
                seq_name = trimmed[1..].to_string();
                state = ParseState::Seed;
            }
            ParseState::Seed => {
                let tokens: Vec<&str> = trimmed.splitn(5, ':').collect();
                let seed = Seed {
                    seq_start: tokens[0].parse::<usize>()?,
                    seq_end: tokens[1].parse::<usize>()?,
                    prf_start: tokens[2].parse::<usize>()?,
                    prf_end: tokens[3].parse::<usize>()?,
                    score: tokens[4].parse::<f32>()?,
                };
                prf_seeds.insert(seq_name.clone(), seed);
                state = ParseState::SeqName;
            }
            ParseState::End => {
                seeds.insert(prf_name.clone(), prf_seeds);
                prf_seeds = HashMap::new();
                state = ParseState::PrfName;
            }
        }

        line.clear();
    }

    Ok(seeds)
}

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

#[derive(Clone)]
pub struct MaxSeedStage {
    seeds: SeedMap,
}

impl MaxSeedStage {
    pub fn new(queries: &Queries, targets: &Fasta) -> Self {
        let queries_and_lengths: Vec<(String, usize)> = match queries {
            Queries::Sequence(fasta) => fasta.iter().map(|s| (s.name.clone(), s.length)).collect(),
            Queries::Profile(p7hmm) => p7hmm.iter().map(|p| (p.name.clone(), p.length)).collect(),
        };

        let seeds: SeedMap = queries_and_lengths
            .into_iter()
            .map(|(q, l)| {
                (
                    q.clone(),
                    targets
                        .iter()
                        .map(|t| {
                            (
                                t.name.clone(),
                                Seed {
                                    seq_start: 1,
                                    seq_end: t.length,
                                    prf_start: 1,
                                    prf_end: l,
                                    score: 0.0,
                                },
                            )
                        })
                        .collect(),
                )
            })
            .collect();

        Self { seeds }
    }
}

impl SeedStage for MaxSeedStage {
    fn run(&mut self, profile: &Profile) -> Option<&HashMap<String, Seed>> {
        self.seeds.get(&profile.name)
    }
}
