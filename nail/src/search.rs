use std::cell::RefCell;
use std::collections::HashMap;
use std::fs::File;
use std::io::{stdout, BufReader};
use std::path::Path;
use std::sync::{Arc, Mutex};
use std::time::{Duration, Instant};

use crate::args::{SearchArgs, SeedMode, SeedsFormat};
use crate::io::{parse_position, write_named_seeds, Seeds};
use crate::io::{Fasta, P7Hmm};
use crate::mmseqs::MmseqsDbPaths;
use crate::pipeline::{
    positional_seeds_path, seed_progressive, seed_static, AlignmentOutput, DefaultAlignStage,
    DefaultCloudSearchStage, FullDpCloudSearchStage, OutputStage, Pipeline, TableOutput,
    BLAST_COLUMNS, DEFAULT_COLUMNS,
};
use crate::stats::{SerialTimed, SetupTimed, Stats, ThreadedTimed};
use crate::util::{guess_query_format_from_query_file, FileFormat};
use crate::util::{term::*, PathExt};

use anyhow::{bail, Context};
use libnail::output::output_tabular::TableFormat;
use libnail::structs::Profile;
use rayon::iter::ParallelIterator;
use rayon::slice::ParallelSlice;
use thread_local::ThreadLocal;

pub enum Queries {
    Sequence(Fasta),
    Profile(Arc<Vec<Profile>>),
}

impl Queries {
    pub fn len(&self) -> usize {
        match self {
            Queries::Sequence(q) => q.len(),
            Queries::Profile(q) => q.len(),
        }
    }

    pub fn name(&self, n: usize) -> anyhow::Result<String> {
        match self {
            Queries::Sequence(q) => q.name(n),
            Queries::Profile(q) => Ok(q[n].name.clone()),
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
            let p7hmm = P7Hmm::from_path(&path).context("failed to index query hmm")?;
            let profiles = p7hmm
                .par_iter()
                .collect::<anyhow::Result<Vec<_>>>()
                .context("failed to parse query hmm")?;
            Ok(Queries::Profile(Arc::new(profiles)))
        }
        _ => {
            panic!()
        }
    }
}

/// Read the seeds file given to --seeds, in the format --seeds-format names.
fn read_seeds(
    path: &Path,
    queries: &Queries,
    targets: &Fasta,
    args: &SearchArgs,
) -> anyhow::Result<Seeds> {
    let max_seeds = args.seed_args.max_seeds;

    match args.io_args.seeds_format {
        SeedsFormat::Positions => {
            Seeds::from_path(path, max_seeds, &parse_position, &parse_position)
        }
        SeedsFormat::Names => {
            let query_positions = (0..queries.len())
                .map(|n| Ok((queries.name(n)?, n)))
                .collect::<anyhow::Result<HashMap<String, usize>>>()?;

            let query_position = |name: &str| {
                query_positions
                    .get(name)
                    .copied()
                    .with_context(|| format!("seed field \"{name}\" is not a query name"))
            };
            let target_position = |name: &str| {
                targets
                    .position(name)?
                    .with_context(|| format!("seed field \"{name}\" is not a target name"))
            };

            Seeds::from_path(path, max_seeds, &query_position, &target_position)
        }
    }
}

pub fn seed(
    queries: &Queries,
    targets: &Fasta,
    stats: &mut Stats,
    args: &mut SearchArgs,
) -> anyhow::Result<Seeds> {
    match args.io_args.seeds_input_path {
        Some(ref path) => read_seeds(path, queries, targets, args),
        None => {
            let now = Instant::now();

            let db_paths = MmseqsDbPaths::new(
                args.io_args
                    .tmp_dir_path
                    .as_ref()
                    .context("args.io_args.tmp_dir_path is somehow unset")?,
            );

            if args.io_args.allow_overwrite {
                db_paths
                    .destroy()
                    .context("failed to remove existing mmseqs DBs")?;
            } else {
                db_paths.check().context("mmseqs DB check failed")?;
            }

            let seeds = match args.seed_args.seed_mode {
                SeedMode::Static => {
                    seed_static(queries, targets, &db_paths, stats, args).context("seeding failed")
                }
                SeedMode::Prog => seed_progressive(queries, targets, &db_paths, stats, args)
                    .context("progessive seeding failed"),
            }?;

            if let (Some(path), SeedsFormat::Names) =
                (&args.io_args.seeds_output_path, args.io_args.seeds_format)
            {
                let positional = BufReader::new(File::open(positional_seeds_path(args)?)?);
                write_named_seeds(
                    positional,
                    path.open(args.io_args.allow_overwrite)?,
                    &|n| queries.name(n),
                    &|n| targets.name(n),
                )
                .context("failed to write named seeds")?;
            }

            stats.set_serial_time(SerialTimed::Seeding, now.elapsed());

            let mut counts: HashMap<String, u64> = HashMap::new();
            for seed in &seeds.seeds {
                *counts.entry(queries.name(seed.prf_pos)?).or_insert(0) += 1;
            }

            stats.set_seed_counts(counts);

            Ok(seeds)
        }
    }
}

pub fn build_pipeline(
    queries: &Queries,
    targets: Fasta,
    stats: Stats,
    args: &mut SearchArgs,
) -> anyhow::Result<Pipeline> {
    let profiles: Arc<Vec<Profile>> = match queries {
        Queries::Sequence(fasta) => Arc::new(
            fasta
                .par_iter()
                .map(|s| Profile::from_blosum_62_and_seq(&s?))
                .collect::<anyhow::Result<Vec<_>>>()
                .context("failed to build profiles")?,
        ),
        Queries::Profile(profiles) => Arc::clone(profiles),
    };

    let mut output = OutputStage::new(args.pipeline_args.e_value_threshold)?;

    if args.ali_to_stdout {
        output.add(AlignmentOutput::new(stdout()))
    } else if let Some(path) = &args.io_args.ali_results_path {
        output.add(AlignmentOutput::new(path.open(true)?))
    }

    if let Some(path) = &args.io_args.tbl_results_path {
        let format = match args.io_args.tbl_format {
            crate::args::TableFormat::Nail => TableFormat::new(&DEFAULT_COLUMNS),
            crate::args::TableFormat::Blast => TableFormat::new(&BLAST_COLUMNS),
        }
        .context("failed to build TableFormat")?;

        output.add(TableOutput::new(path.open(true)?, format));
    }

    Ok(Pipeline {
        profiles,
        prf: None,
        targets: Arc::new(targets),
        cloud_search: match args.dev_args.full_dp {
            true => Box::<FullDpCloudSearchStage>::default(),
            false => Box::new(DefaultCloudSearchStage::new(args)),
        },
        align: Box::new(
            DefaultAlignStage::new(args).context("failed to create DefaultAlignStage")?,
        ),
        output: Arc::new(Mutex::new(output)),
        stats,
    })
}

fn timed_step<T>(
    label: &str,
    step: impl FnOnce() -> anyhow::Result<T>,
) -> anyhow::Result<(T, Duration)> {
    let now = Instant::now();
    println!("{label}...");
    let result = step()?;
    let elapsed = now.elapsed();
    println!(
        "\x1b[A{:<28}done ({:.2}s)",
        format!("{label}..."),
        elapsed.as_secs_f64()
    );
    Ok((result, elapsed))
}

pub fn search(mut args: SearchArgs) -> anyhow::Result<()> {
    let start_time = Instant::now();

    let (queries, query_index_time) =
        timed_step("indexing query database", || read_queries(&args.query_path))?;

    let (targets, target_index_time) = timed_step("indexing target database", || {
        Fasta::from_path(&args.target_path).context("failed to index target fasta")
    })?;

    let mut stats = Stats::new(&queries, targets.len())?;
    stats.set_num_threads(args.num_threads);
    stats.set_setup_time(SetupTimed::QueryIndex, query_index_time);
    stats.set_setup_time(SetupTimed::TargetIndex, target_index_time);

    match args.expert_args.target_database_size {
        Some(_) => {}
        None => args.expert_args.target_database_size = Some(targets.len()),
    }

    if args.pipeline_args.only_seed {
        timed_step("seeding", || {
            seed(&queries, &targets, &mut stats, &mut args)
        })?;
        return Ok(());
    }

    let (mut pipeline, build_time) = timed_step("building pipeline", || {
        build_pipeline(&queries, targets, stats, &mut args).context("failed to build pipeline")
    })?;
    pipeline
        .stats
        .set_setup_time(SetupTimed::PipelineBuild, build_time);
    pipeline
        .stats
        .set_setup_time(SetupTimed::Total, start_time.elapsed());
    pipeline
        .stats
        .set_serial_time(SerialTimed::Setup, start_time.elapsed());

    let (seeds, _) = timed_step("seeding", || {
        seed(&queries, &pipeline.targets, &mut pipeline.stats, &mut args)
    })?;

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

            let lock_now = Instant::now();
            let (output_stats, lock_time) = match pipeline.output.lock() {
                Ok(mut guard) => {
                    let lock_time = lock_now.elapsed();
                    (guard.run(&res)?, lock_time)
                }
                Err(_) => bail!("mutex poisoned"),
            };

            pipeline
                .stats
                .add_threaded_time(ThreadedTimed::OutputMutex, lock_time);

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
        "\x1b[A{:<28}done ({:.2}s)\n",
        "running nail pipeline...",
        align_timer.elapsed().as_secs_f64()
    );

    if let Some(path) = &args.io_args.tbl_results_path {
        println!("tabular results written to: {BLUE}{path:?}{RESET}");
    }

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
