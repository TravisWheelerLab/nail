use std::{
    collections::HashMap,
    fs::File,
    io::{BufWriter, Write},
    path::Path,
};

use libnail::{align::structs::Seed, structs::Profile};

use crate::{
    args::SearchArgs,
    io::{Database, Fasta, P7Hmm, SeedList, Seeds},
    mmseqs::{
        run_mmseqs_search, write_mmseqs_profile_database, write_mmseqs_sequence_database,
        MmseqsDbPaths, MmseqsScoreModel,
    },
};

fn merge_seeds<P: AsRef<Path>>(
    seeds_a: Seeds,
    mut seeds_b: Seeds,
    path: P,
) -> anyhow::Result<Seeds> {
    let mut out = BufWriter::new(File::create(path.as_ref())?);
    seeds_a
        .iter()
        .map(|(prf, a)| (prf, a, seeds_b.get(prf).unwrap_or_default()))
        .try_for_each(|(prf, a, b)| {
            let mut hash: HashMap<String, Seed> = HashMap::new();
            a.into_iter()
                .chain(b)
                .for_each(|(seq, seed)| match hash.get(&seq) {
                    Some(existing) => {
                        if existing.score < seed.score {
                            hash.insert(seq, seed);
                        }
                    }
                    None => {
                        hash.insert(seq, seed);
                    }
                });

            hash.into_iter().try_for_each(|(seq, seed)| {
                writeln!(
                    out,
                    "{prf} {seq} {} {} {} {} {}",
                    seed.prf_start, seed.prf_end, seed.seq_start, seed.seq_end, seed.score
                )
            })
        })?;

    // we need to explicitly drop the file
    // handle before reading from it
    drop(out);

    Seeds::from_path(path)
}

pub fn seed_profile_to_sequence(
    profiles: &P7Hmm,
    seqs: &Fasta,
    args: &SearchArgs,
) -> anyhow::Result<Seeds> {
    let paths = MmseqsDbPaths::new(&args.io_args.temp_dir_path);

    write_mmseqs_sequence_database(seqs, &paths.target_db)?;

    let align_tsv_a = &paths.dir()?.join("align_a.tsv");
    write_mmseqs_profile_database(profiles, &paths.query_db, None)?;
    run_mmseqs_search(&paths, align_tsv_a, args, MmseqsScoreModel::Profile)?;
    let seeds_a = Seeds::from_path(align_tsv_a)?;

    if !args.pipeline_args.double_seed {
        return Ok(seeds_a);
    }

    let align_tsv_b = &paths.dir()?.join("align_b.tsv");
    write_mmseqs_profile_database(profiles, &paths.query_db, Some(1.0))?;
    run_mmseqs_search(&paths, align_tsv_b, args, MmseqsScoreModel::Profile)?;
    let seeds_b = Seeds::from_path(align_tsv_b)?;

    let align_tsv_merge = &paths.dir()?.join("align_merge.tsv");
    merge_seeds(seeds_a, seeds_b, align_tsv_merge)
}

pub fn seed_sequence_to_sequence(
    queries: &Fasta,
    targets: &Fasta,
    args: &SearchArgs,
) -> anyhow::Result<Seeds> {
    let paths = MmseqsDbPaths::new(&args.io_args.temp_dir_path);

    write_mmseqs_sequence_database(targets, &paths.target_db)?;
    write_mmseqs_sequence_database(queries, &paths.query_db)?;

    let align_tsv_a = &paths.dir()?.join("align_a.tsv");
    run_mmseqs_search(&paths, align_tsv_a, args, MmseqsScoreModel::Blosum62)?;
    let seeds_a = Seeds::from_path(align_tsv_a)?;

    if !args.pipeline_args.double_seed {
        return Ok(seeds_a);
    }

    let align_tsv_b = &paths.dir()?.join("align_b.tsv");
    run_mmseqs_search(&paths, align_tsv_b, args, MmseqsScoreModel::Blosum80)?;
    let seeds_b = Seeds::from_path(align_tsv_b)?;

    let align_tsv_merge = &paths.dir()?.join("align_merge.tsv");
    merge_seeds(seeds_a, seeds_b, align_tsv_merge)
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
