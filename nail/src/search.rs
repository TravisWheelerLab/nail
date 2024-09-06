use std::collections::HashMap;
use std::fs::File;
use std::io::{stdout, BufReader, BufWriter};
use std::path::Path;
use std::sync::{Arc, Mutex};
use std::time::Instant;

use crate::args::{SearchArgs, SeedArgs};
use crate::pipeline::{
    run_pipeline_profile_to_sequence, run_pipeline_sequence_to_sequence, seed_profile_to_sequence,
    seed_sequence_to_sequence, DefaultAlignStep, DefaultCloudSearchStep, DefaultSeedStep,
    FullDpCloudSearchStep, OutputStep, Pipeline, SeedMap,
};
use crate::stats::Stats;
use crate::util::{guess_query_format_from_query_file, FileFormat, PathBufExt};

use libnail::structs::{Hmm, Profile, Sequence};

use anyhow::Context;
use serde::Serialize;

pub enum Queries {
    Sequence(Vec<Sequence>),
    Profile(Vec<Profile>),
}

fn read_queries(path: impl AsRef<Path>) -> anyhow::Result<Queries> {
    let query_format = guess_query_format_from_query_file(&path)?;

    match query_format {
        FileFormat::Fasta => {
            let queries =
                Sequence::amino_from_fasta(&path).context("failed to read query fasta")?;

            Ok(Queries::Sequence(queries))
        }
        FileFormat::Hmm => {
            let hmm_file = File::open(&path)?;
            let queries: Vec<Profile> = Hmm::from_p7hmm(hmm_file)
                .context("failed to read query hmm")?
                .iter()
                .map(Profile::new)
                .collect();

            Ok(Queries::Profile(queries))
        }
        _ => {
            panic!()
        }
    }
}

pub fn seed(args: SeedArgs) -> anyhow::Result<()> {
    let queries = read_queries(&args.query_path)?;
    let targets =
        Sequence::amino_from_fasta(&args.target_path).context("failed to read target fasta")?;

    let seeds = match queries {
        Queries::Sequence(ref queries) => seed_sequence_to_sequence(
            queries,
            &targets,
            args.common_args.num_threads,
            &args.mmseqs_args,
        )?,
        Queries::Profile(ref queries) => seed_profile_to_sequence(
            queries,
            &targets,
            args.common_args.num_threads,
            &args.mmseqs_args,
        )?,
    };

    let writer = BufWriter::new(args.seeds_path.open(true)?);
    let mut serializer = serde_json::Serializer::new(writer);
    seeds.serialize(&mut serializer)?;

    Ok(())
}

pub fn search(mut args: SearchArgs) -> anyhow::Result<()> {
    {
        // quickly make sure we can write the results
        args.output_args.tbl_results_path.open(true)?;
        if let Some(path) = &args.output_args.ali_results_path {
            path.open(true)?;
        }
    }

    let mut stats = Stats::default();

    let now = Instant::now();
    let queries = read_queries(&args.query_path)?;
    stats.add_input_time(now.elapsed());

    let now = Instant::now();
    let targets =
        Sequence::amino_from_fasta(&args.target_path).context("failed to read target fasta")?;
    stats.add_input_time(now.elapsed());

    match args.nail_args.target_database_size {
        Some(_) => {}
        None => args.nail_args.target_database_size = Some(targets.len()),
    }

    let seeds = match args.seeds_path {
        Some(ref path) => {
            let mut seeds: SeedMap = HashMap::new();

            let reader = BufReader::new(std::fs::File::open(path)?);
            let stream = serde_json::Deserializer::from_reader(reader);

            let now = Instant::now();
            for entry in stream.into_iter::<SeedMap>() {
                let entry = entry?;
                seeds.extend(entry);
            }
            stats.add_input_time(now.elapsed());

            seeds
        }
        None => {
            let now = Instant::now();
            let seeds = match queries {
                Queries::Sequence(ref queries) => seed_sequence_to_sequence(
                    queries,
                    &targets,
                    args.common_args.num_threads,
                    &args.mmseqs_args,
                )?,
                Queries::Profile(ref queries) => seed_profile_to_sequence(
                    queries,
                    &targets,
                    args.common_args.num_threads,
                    &args.mmseqs_args,
                )?,
            };
            stats.set_mmseqs_time(now.elapsed());

            seeds
        }
    };

    let mut pipeline = Pipeline {
        seed: Box::new(DefaultSeedStep::new(seeds)),
        cloud_search: match args.nail_args.full_dp {
            true => Box::<FullDpCloudSearchStep>::default(),
            false => Box::new(DefaultCloudSearchStep::new(&args)),
        },
        align: Box::new(DefaultAlignStep::new(&args)),
        output: Arc::new(Mutex::new(OutputStep::new(&args.output_args)?)),
        stats,
    };

    let target_map: HashMap<String, Sequence> =
        targets.into_iter().map(|t| (t.name.clone(), t)).collect();

    match queries {
        Queries::Sequence(queries) => {
            run_pipeline_sequence_to_sequence(&queries, &target_map, pipeline);
        }
        Queries::Profile(mut queries) => {
            run_pipeline_profile_to_sequence(&mut queries, &target_map, &mut pipeline);
        }
    }

    Ok(())
}
