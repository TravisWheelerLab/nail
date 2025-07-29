use std::{
    collections::HashMap,
    io::{BufRead, BufReader, BufWriter, Read, Seek, SeekFrom, Write},
    sync::Arc,
};

use libnail::{align::structs::Seed, structs::Profile};

use crate::{
    args::SearchArgs,
    io::{ByteBufferExt, Fasta, P7Hmm, ProfileDatabase, ReadSeekExt, ReadState, SequenceDatabase},
    mmseqs::{
        run_mmseqs_search, seeds_from_mmseqs_align_tsv, write_mmseqs_profile_database,
        write_mmseqs_sequence_database, MmseqsDbPaths, MmseqsScoreModel,
    },
    search::Queries,
};

use anyhow::anyhow;

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

pub fn read_seed_map<R: Read + Seek>(data: &mut R) -> anyhow::Result<SeedMap> {
    let mut seeds: SeedMap = HashMap::new();

    #[derive(Debug)]
    enum ParseState {
        PrfName,
        SeqName,
        Coords,
        End,
    }
    use ParseState::*;

    // 2<<16 == 2^17 | ~128KiB buffer
    let mut buf = vec![0u8; 2 << 16];

    // search for the first newline in the first 2^13 bytes
    let first_newline_pos = data.read_to_first_newline(&mut buf[0..2 << 12])?;

    let mut state = SeqName;
    let mut prf_name = (&buf[0..first_newline_pos])
        .str(0, first_newline_pos - 1)?
        .to_string();
    let mut seq_name = "".to_string();
    let mut prf_seeds = HashMap::new();
    data.seek(SeekFrom::Start(first_newline_pos as u64))?;

    while let Ok(read_state) = data.read_with_state(&mut buf) {
        let buf_slice = match read_state {
            ReadState::Reading(n) => {
                let last_newline_pos = buf[..n]
                    .iter()
                    .rposition(|&b| b == b'\n')
                    .ok_or(anyhow!(""))?;

                let n_bytes_back = (n - last_newline_pos) as i64;
                data.seek_relative(-n_bytes_back)?;
                &buf[..last_newline_pos]
            }
            ReadState::Final(n) => &buf[..n],
            ReadState::Done => break,
        };

        let mut i = 0;
        while i < buf_slice.len() {
            let word = buf_slice.word_from(i + 1)?;

            i += word.len() + 1;

            match state {
                PrfName => {
                    prf_name.push_str(word);
                    state = SeqName;
                }
                SeqName => {
                    if word == "//" {
                        state = End;
                        i -= 1;
                        continue;
                    }
                    seq_name.push_str(&word[1..]);
                    state = Coords;
                }
                Coords => {
                    let mut tokens = word.split(':');
                    let seed = Seed {
                        seq_start: tokens.next().unwrap().parse::<usize>()?,
                        seq_end: tokens.next().unwrap().parse::<usize>()?,
                        prf_start: tokens.next().unwrap().parse::<usize>()?,
                        prf_end: tokens.next().unwrap().parse::<usize>()?,
                        score: tokens.next().unwrap().parse::<f32>()?,
                    };
                    seq_name.shrink_to_fit();
                    prf_seeds.insert(seq_name.clone(), seed);
                    seq_name.clear();
                    state = SeqName;
                }
                End => {
                    prf_name.shrink_to_fit();
                    seeds.insert(prf_name.clone(), prf_seeds);
                    prf_name.clear();
                    prf_seeds = HashMap::new();
                    state = PrfName;
                }
            }
        }
    }

    Ok(seeds)
}

#[derive(Default, Clone)]
pub struct DefaultSeedStage {
    seeds: Arc<SeedMap>,
}

impl DefaultSeedStage {
    pub fn new(seeds: SeedMap) -> Self {
        DefaultSeedStage {
            seeds: Arc::new(seeds),
        }
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
