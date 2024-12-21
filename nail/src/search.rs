use std::collections::HashMap;
use std::fs::File;
use std::io::{stdout, BufReader, BufWriter};
use std::path::Path;
use std::time::Instant;

use crate::args::{SearchArgs, SeedArgs};
use crate::io::Fasta;
use crate::pipeline::{
    run_pipeline_profile_to_sequence, run_pipeline_sequence_to_sequence, seed_profile_to_sequence,
    seed_sequence_to_sequence, DefaultAlignStage, DefaultCloudSearchStage, DefaultSeedStage,
    FullDpCloudSearchStage, OutputStage, Pipeline, SeedMap,
};
use crate::stats::{SerialTimed, Stats};
use crate::util::{guess_query_format_from_query_file, FileFormat, PathBufExt};

use libnail::structs::{Hmm, Profile};

use anyhow::Context;
use serde::Serialize;

pub enum Queries {
    Sequence(Fasta),
    Profile(Vec<Profile>),
}

impl Queries {
    pub fn len(&self) -> usize {
        match self {
            Queries::Sequence(q) => q.len(),
            Queries::Profile(q) => q.len(),
        }
    }
}

fn read_queries(path: impl AsRef<Path>) -> anyhow::Result<Queries> {
    let query_format = guess_query_format_from_query_file(&path)?;

    match query_format {
        FileFormat::Fasta => {
            let queries = Fasta::from_path(&path).context("failed to read query fasta")?;
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
    let targets = Fasta::from_path(&args.target_path).context("failed to read target fasta")?;

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
    let start_time = Instant::now();
    {
        // quickly make sure we can write the results
        args.output_args.tbl_results_path.open(true)?;
        if let Some(path) = &args.output_args.ali_results_path {
            path.open(true)?;
        }
    }

    let now = Instant::now();
    println!("reading query database...");
    let queries = read_queries(&args.query_path)?;
    println!(
        "\x1b[Areading query database...   done ({:.2}s)",
        now.elapsed().as_secs_f64()
    );

    let now = Instant::now();
    println!("indexing target database...");
    let targets = Fasta::from_path(&args.target_path).context("failed to read target fasta")?;
    println!(
        "\x1b[Aindexing target database... done ({:.2}s)",
        now.elapsed().as_secs_f64()
    );

    let mut stats = Stats::new(&queries, &targets);

    match args.nail_args.target_database_size {
        Some(_) => {}
        None => args.nail_args.target_database_size = Some(targets.len()),
    }

    let seeds = match args.seeds_path {
        Some(ref path) => {
            let mut seeds: SeedMap = HashMap::new();

            let reader = BufReader::new(std::fs::File::open(path)?);
            let stream = serde_json::Deserializer::from_reader(reader);

            for entry in stream.into_iter::<SeedMap>() {
                let entry = entry?;
                seeds.extend(entry);
            }

            seeds
        }
        None => {
            let now = Instant::now();
            println!("running mmseqs...");
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
            stats.set_serial_time(SerialTimed::Seeding, now.elapsed());
            println!(
                "\x1b[Arunning mmseqs...           done ({:.2}s)",
                now.elapsed().as_secs_f64()
            );
            seeds
        }
    };

    let mut pipeline = Pipeline {
        targets,
        seed: Box::new(DefaultSeedStage::new(seeds)),
        cloud_search: match args.nail_args.full_dp {
            true => Box::<FullDpCloudSearchStage>::default(),
            false => Box::new(DefaultCloudSearchStage::new(&args)),
        },
        align: Box::new(DefaultAlignStage::new(&args)),
        output: OutputStage::new(&args.output_args)?,
        stats,
    };

    println!("running nail pipeline...");
    let align_timer = Instant::now();
    match queries {
        Queries::Sequence(queries) => {
            run_pipeline_sequence_to_sequence(&queries, &mut pipeline);
        }
        Queries::Profile(mut queries) => {
            run_pipeline_profile_to_sequence(&mut queries, &mut pipeline);
        }
    }

    pipeline
        .stats
        .set_serial_time(SerialTimed::Alignment, align_timer.elapsed());

    println!(
        "\x1b[Arunning nail pipeline...    done ({:.2}s)\n",
        align_timer.elapsed().as_secs_f64()
    );

    pipeline
        .stats
        .set_serial_time(SerialTimed::Total, start_time.elapsed());

    // pipeline.stats.write(&mut stdout())?;

    Ok(())
}
