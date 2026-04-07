use std::cell::RefCell;
use std::collections::HashMap;
use std::io::stdout;
use std::path::Path;
use std::sync::Arc;
use std::time::Instant;

use crate::args::SearchArgs;
use crate::io::{Fasta, P7Hmm, Seeds2};
use crate::mmseqs::MmseqsDbPaths;
use crate::pipeline::{
    seed_max_seqs, seed_progressive, DefaultAlignStage, DefaultCloudSearchStage,
    FullDpCloudSearchStage, OutputStage, Pipeline,
};
use crate::stats::{SerialTimed, Stats, ThreadedTimed};
use crate::util::{guess_query_format_from_query_file, FileFormat};

use anyhow::Context;
use libnail::structs::Profile;
use rayon::iter::{IntoParallelRefIterator, ParallelIterator};
use rayon::slice::ParallelSlice;
use thread_local::ThreadLocal;

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

pub fn read_queries(path: impl AsRef<Path>) -> anyhow::Result<Queries> {
    let query_format = guess_query_format_from_query_file(&path)?;

    match query_format {
        FileFormat::Fasta => {
            let queries = Fasta::from_path(&path).context("failed to index query fasta")?;
            Ok(Queries::Sequence(queries))
        }
        FileFormat::Hmm => {
            let queries = P7Hmm::from_path(&path).context("failed to index query hmm")?;
            Ok(Queries::Profile(queries))
        }
        _ => {
            panic!()
        }
    }
}

pub fn seed(
    queries: &Queries,
    targets: &Fasta,
    stats: &mut Stats,
    args: &mut SearchArgs,
) -> anyhow::Result<Seeds2> {
    match args.io_args.seeds_input_path {
        Some(ref path) => Seeds2::from_path(path),
        None => {
            let now = Instant::now();

            let db_paths = MmseqsDbPaths::new(&args.io_args.temp_dir_path);
            if args.io_args.allow_overwrite {
                db_paths
                    .destroy()
                    .context("failed to remove existing mmseqs DBs")?;
            } else {
                db_paths.check().context("mmseqs DB check failed")?;
            }

            let seeds = if args.mmseqs_args.prog_seed {
                seed_progressive(queries, targets, &db_paths, stats, args)
                    .context("progessive seeding failed")
            } else {
                seed_max_seqs(queries, targets, &db_paths, stats, args).context("seeding failed")
            }?;

            stats.set_serial_time(SerialTimed::Seeding, now.elapsed());

            let mut counts: HashMap<String, u64> = HashMap::new();
            for seed in &seeds.seeds {
                *counts.entry(seed.prf.clone()).or_insert(0) += 1;
            }

            stats.set_seed_counts(counts);

            Ok(seeds)
        }
    }
}

pub fn build_pipeline(
    queries: Queries,
    targets: Fasta,
    stats: Stats,
    args: &mut SearchArgs,
) -> anyhow::Result<Pipeline> {
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

    Ok(Pipeline {
        profiles,
        prf: None,
        targets,
        cloud_search: match args.dev_args.full_dp {
            true => Box::<FullDpCloudSearchStage>::default(),
            false => Box::new(DefaultCloudSearchStage::new(args)),
        },
        align: Box::new(
            DefaultAlignStage::new(args).context("failed to create DefaultAlignStage")?,
        ),
        output: OutputStage::new(args).context("failed to create OutputStage")?,
        stats,
    })
}

pub fn search(mut args: SearchArgs) -> anyhow::Result<()> {
    let start_time = Instant::now();

    let now = Instant::now();
    println!("indexing query database...");
    let queries = read_queries(&args.query_path)?;
    println!(
        "\x1b[Aindexing query database...  done ({:.2}s)",
        now.elapsed().as_secs_f64()
    );

    let now = Instant::now();
    println!("indexing target database...");
    let targets = Fasta::from_path(&args.target_path).context("failed to index target fasta")?;
    println!(
        "\x1b[Aindexing target database... done ({:.2}s)",
        now.elapsed().as_secs_f64()
    );

    let mut stats = Stats::new(queries.len(), targets.len());

    match args.expert_args.target_database_size {
        Some(_) => {}
        None => args.expert_args.target_database_size = Some(targets.len()),
    }

    let now = Instant::now();
    println!("seeding...");
    let seeds = seed(&queries, &targets, &mut stats, &mut args)?;
    println!(
        "\x1b[Aseeding...                  done ({:.2}s)",
        now.elapsed().as_secs_f64()
    );

    if args.pipeline_args.only_seed {
        return Ok(());
    }

    let mut pipeline = build_pipeline(queries, targets, stats, &mut args)?;

    let align_timer = Instant::now();
    println!("running nail pipeline...");

    let tl_pipeline: ThreadLocal<RefCell<Pipeline>> = ThreadLocal::new();

    seeds
        .seeds
        .par_chunks(100)
        .panic_fuse()
        .try_for_each(|chunk| -> anyhow::Result<()> {
            let now = Instant::now();
            let mut pipeline = tl_pipeline
                .get_or(|| RefCell::new(pipeline.clone()))
                .borrow_mut();

            let res = chunk
                .iter()
                .map(|seed| pipeline.run(seed))
                .collect::<Result<Vec<_>, _>>()?;

            let output_stats = pipeline.output.run(&res)?;
            pipeline.stats.add_sample(&res, &output_stats);

            pipeline
                .stats
                .add_threaded_time(ThreadedTimed::Total, now.elapsed());

            Ok(())
        })?;

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

    pipeline.stats.write_max_seqs_report(&args)?;

    Ok(())
}
