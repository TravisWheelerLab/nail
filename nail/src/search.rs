use std::collections::HashMap;
use std::io::stdout;
use std::path::{Path, PathBuf};
use std::str::FromStr;
use std::sync::Arc;
use std::time::Instant;

use crate::args::SearchArgs;
use crate::io::{Fasta, P7Hmm, Seeds2};
use crate::pipeline::{
    run_pipeline_profile_to_sequence, seed_profile_to_sequence,
    seed_profile_to_sequence_progressive, seed_progressive, seed_sequence_to_sequence,
    DefaultAlignStage, DefaultCloudSearchStage, FullDpCloudSearchStage, OutputStage, Pipeline,
};
use crate::stats::{SerialTimed, Stats};
use crate::util::{guess_query_format_from_query_file, FileFormat, PathBufExt};

use anyhow::Context;
use libnail::structs::Profile;
use rayon::iter::ParallelIterator;

pub enum Queries {
    Sequence(Fasta),
    Profile(P7Hmm),
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
            let queries = P7Hmm::from_path(&path).context("failed to open query hmm")?;
            Ok(Queries::Profile(queries))
        }
        _ => {
            panic!()
        }
    }
}

pub fn search(mut args: SearchArgs) -> anyhow::Result<()> {
    let start_time = Instant::now();

    if args.ali_to_stdout {
        args.io_args.tbl_results_path = None
    }

    if args.pipeline_args.only_seed && args.io_args.seeds_output_path.is_none() {
        args.io_args.seeds_output_path = Some(PathBuf::from_str("./seeds.json")?);
    }

    {
        // quickly make sure we can write to all of the results paths
        if let Some(path) = &args.io_args.tbl_results_path {
            path.open(args.io_args.allow_overwrite)?;
        }

        if let Some(path) = &args.io_args.ali_results_path {
            path.open(args.io_args.allow_overwrite)?;
        }

        if let Some(path) = &args.io_args.seeds_output_path {
            path.open(args.io_args.allow_overwrite)?;
        }

        if let Some(path) = &args.dev_args.stats_results_path {
            path.open(args.io_args.allow_overwrite)?;
        }
    }

    let now = Instant::now();
    println!("indexing query database...");
    let queries = read_queries(&args.query_path)?;
    println!(
        "\x1b[Aindexing query database...  done ({:.2}s)",
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

    match args.expert_args.target_database_size {
        Some(_) => {}
        None => args.expert_args.target_database_size = Some(targets.len()),
    }

    let seeds = match args.io_args.seeds_input_path {
        Some(ref path) => {
            let now = Instant::now();
            println!("reading seeds...");
            let seeds = Seeds2::from_path(path);
            println!(
                "\x1b[Areading seeds...            done ({:.2}s)",
                now.elapsed().as_secs_f64()
            );

            seeds
        }
        None => {
            let now = Instant::now();
            println!("running mmseqs...");
            let seeds = seed_progressive(&queries, &targets, &mut stats, &args);
            stats.set_serial_time(SerialTimed::Seeding, now.elapsed());
            println!(
                "\x1b[Arunning mmseqs...           done ({:.2}s)",
                now.elapsed().as_secs_f64()
            );

            seeds
        }
    }?;

    if args.pipeline_args.only_seed {
        return Ok(());
    }

    let profiles: Arc<HashMap<String, Profile>> = Arc::new(
        match queries {
            Queries::Sequence(fasta) => fasta
                .par_iter()
                .filter_map(|s| Profile::from_blosum_62_and_seq(&s).ok())
                .collect::<Vec<_>>(),
            Queries::Profile(p7hmm) => p7hmm.par_iter().collect::<Vec<_>>(),
        }
        .into_iter()
        .map(|p| (p.name.clone(), p))
        .collect(),
    );

    let mut pipeline = Pipeline {
        profiles,
        prf: None,
        targets,
        cloud_search: match args.dev_args.full_dp {
            true => Box::<FullDpCloudSearchStage>::default(),
            false => Box::new(DefaultCloudSearchStage::new(&args)),
        },
        align: Box::new(
            DefaultAlignStage::new(&args).context("failed to create DefaultAlignStage")?,
        ),
        output: OutputStage::new(&args).context("failed to create OutputStage")?,
        stats,
    };

    let align_timer = Instant::now();
    println!("running nail pipeline...");

    run_pipeline_profile_to_sequence(&mut pipeline, seeds);

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

    if args.print_summary_stats {
        args.write(&mut stdout())?;
        pipeline.stats.write(&mut stdout())?;
    }

    Ok(())
}
