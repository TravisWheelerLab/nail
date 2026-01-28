use std::{io::Write, time::Instant};

use anyhow::Context;
use libnail::align::{structs::Seed, Bits};

use crate::{
    args::SearchArgs,
    io::{Database, Fasta, Seeds2},
    mmseqs::{
        consts::{ALIGN_DBTYPE, PREFILTER_DBTYPE},
        run_mmseqs_align, run_mmseqs_convertalis, run_mmseqs_prefilter,
        write_mmseqs_profile_database, write_mmseqs_sequence_database, ByteBuffer, MmseqsDbPaths,
        PrefilterDb,
    },
    pipeline::StageResult,
    search::Queries,
    stats::Stats,
    util::PathBufExt,
};

pub type SeedStageResult = StageResult<Seed, SeedStageStats>;

impl SeedStageResult {}

#[derive(Default)]
pub struct SeedStageStats {
    pub score: Bits,
    pub e_value: f64,
}

pub fn seed_max_seqs(
    queries: &Queries,
    seqs: &Fasta,
    db_paths: &MmseqsDbPaths,
    stats: &mut Stats,
    args: &SearchArgs,
) -> anyhow::Result<Seeds2> {
    let time_start = Instant::now();

    // ---

    write_mmseqs_sequence_database(seqs, &db_paths.target_db)?;
    match queries {
        Queries::Sequence(fasta) => write_mmseqs_sequence_database(fasta, &db_paths.query_db)?,
        Queries::Profile(hmm) => write_mmseqs_profile_database(hmm.values(), &db_paths.query_db)?,
    }

    let now = Instant::now();
    run_mmseqs_prefilter(
        &db_paths.query_db,
        &db_paths.target_db,
        &db_paths.prefilter_db,
        None,
        args,
    )
    .context("mmseqs prefilter failed")?;

    stats.set_mmseqs_time(crate::stats::MmseqsTimed::Prefilter, now.elapsed());

    // ---

    let now = Instant::now();
    run_mmseqs_align(
        &db_paths.query_db,
        &db_paths.target_db,
        &db_paths.prefilter_db,
        &db_paths.align_db,
        None,
        args,
    )
    .context("mmseqs align failed")?;

    stats.add_mmseqs_time(crate::stats::MmseqsTimed::Align, now.elapsed());

    // ---

    let align_tsv = match &args.io_args.seeds_output_path {
        Some(path) => path,
        None => &args.io_args.temp_dir_path.join("seeds.tsv"),
    };

    let now = Instant::now();
    run_mmseqs_convertalis(
        &db_paths.query_db,
        &db_paths.target_db,
        &db_paths.align_db,
        align_tsv,
        args,
    )
    .context("mmseqs convertalis failed")?;

    stats.set_mmseqs_time(crate::stats::MmseqsTimed::Convertalis, now.elapsed());

    // ---

    let now = Instant::now();
    let seeds = Seeds2::from_path(align_tsv).context("failed to build seeds")?;

    stats.set_mmseqs_time(crate::stats::MmseqsTimed::Index, now.elapsed());
    stats.set_mmseqs_time(crate::stats::MmseqsTimed::Total, time_start.elapsed());

    Ok(seeds)
}

pub fn seed_progressive(
    queries: &Queries,
    seqs: &Fasta,
    db_paths: &MmseqsDbPaths,
    stats: &mut Stats,
    args: &SearchArgs,
) -> anyhow::Result<Seeds2> {
    let time_start = Instant::now();

    if !&db_paths
        .prog_dir
        .try_exists()
        .with_context(|| format!("failed to check existence of: {:?}", db_paths.prog_dir))?
    {
        std::fs::create_dir(&db_paths.prog_dir).with_context(|| {
            format!("failed to create prog directory: {:?}", &db_paths.prog_dir)
        })?;
    }

    // ---

    write_mmseqs_sequence_database(seqs, &db_paths.target_db)?;
    match queries {
        Queries::Sequence(fasta) => write_mmseqs_sequence_database(fasta, &db_paths.query_db)?,
        Queries::Profile(hmm) => write_mmseqs_profile_database(hmm.values(), &db_paths.query_db)?,
    }

    let now = Instant::now();
    run_mmseqs_prefilter(
        &db_paths.query_db,
        &db_paths.target_db,
        &db_paths.prefilter_db,
        None,
        args,
    )
    .context("mmseqs prefilter failed")?;

    stats.set_mmseqs_time(crate::stats::MmseqsTimed::Prefilter, now.elapsed());

    // ---

    let mut pdb = PrefilterDb::from_path(&db_paths.prefilter_db)
        .context("failed to open initial prefilter DB")?;

    #[derive(PartialEq, Debug)]
    enum State {
        Active,
        Final(usize),
        Terminated,
    }

    let mut i = 0;
    let mut n_take = args.mmseqs_args.prog_n.context("prog_n unset")?;
    let prog_frac = args.mmseqs_args.prog_f.context("prog_f unset")?;

    let mut state: Vec<State> = (0..queries.len()).map(|_| State::Active).collect();
    let mut prog_adbs = vec![];
    while state.contains(&State::Active) {
        let prog_iter_dir = &db_paths.prog_dir.join(i.to_string());
        let prog_pdb_path = prog_iter_dir.join("pdb");

        std::fs::create_dir(prog_iter_dir).with_context(|| {
            format!("failed to create prog iteration directory: {prog_iter_dir:?}")
        })?;

        {
            // note: scoped to drop file handles and force a write
            let mut prog_pfdb = prog_pdb_path.open(true)?;
            let mut prog_pfdb_index = prog_pdb_path.with_extension("index").open(true)?;
            prog_pdb_path
                .with_extension("dbtype")
                .open(true)?
                .write_all(PREFILTER_DBTYPE)?;

            let mut offset = 0;
            for (prf_idx, prf_state) in state.iter_mut().enumerate() {
                let record_bytes = match prf_state {
                    State::Active => match pdb.next_n(prf_idx, n_take)? {
                        ByteBuffer::Complete(buf) => buf,
                        ByteBuffer::Partial(buf, n_retrieved) => {
                            *prf_state = State::Final(n_retrieved);
                            buf
                        }
                        ByteBuffer::Empty => {
                            *prf_state = State::Terminated;
                            &[]
                        }
                    },
                    State::Terminated => &[],
                    State::Final(_) => unreachable!(),
                };

                prog_pfdb.write_all(record_bytes)?;
                prog_pfdb.write_all(&[0])?;

                let n_written = record_bytes.len() + 1;

                let s = format!("{}\t{offset}\t{}\n", prf_idx, n_written);
                prog_pfdb_index.write_all(s.as_bytes())?;

                offset += n_written;
            }
        }

        let prog_adb_path = prog_iter_dir.join("adb");

        let now = Instant::now();
        run_mmseqs_align(
            &db_paths.query_db,
            &db_paths.target_db,
            &prog_pdb_path,
            &prog_adb_path,
            None,
            args,
        )
        .context("mmseqs align failed")?;

        stats.add_mmseqs_time(crate::stats::MmseqsTimed::Align, now.elapsed());

        let mut prog_adb =
            PrefilterDb::from_path(prog_adb_path).context("failed to open prog align DB")?;

        let mut report_out = prog_pdb_path.with_file_name("report.txt").open(true)?;

        for (prf_idx, prf_state) in state.iter_mut().enumerate() {
            match prf_state {
                State::Active => {
                    let record_bytes = prog_adb.get(prf_idx)?;
                    let mut cnt = 0;
                    for b in record_bytes.iter() {
                        if *b == b'\n' {
                            cnt += 1;
                        }
                    }

                    let frac = cnt as f32 / n_take as f32;

                    if frac < prog_frac {
                        *prf_state = State::Terminated;
                    }
                    writeln!(report_out, "{prf_idx}: {cnt} / {n_take} | {frac:0.3}")?;
                }
                State::Final(n) => {
                    writeln!(report_out, "{prf_idx}: ({n})")?;
                    *prf_state = State::Terminated;
                }
                State::Terminated => continue,
            }
        }

        i += 1;
        n_take *= 2;
        prog_adbs.push(prog_adb);
    }

    // --

    {
        // note: scoped to drop file handles and force a write
        let adb_dir = db_paths
            .align_db
            .parent()
            .context("failed to produce mmseqs align DB directory path")?;

        if !adb_dir
            .try_exists()
            .with_context(|| format!("failed to check existence of: {adb_dir:?}"))?
        {
            std::fs::create_dir(adb_dir).context("failed to create mmseqs align DB directory")?;
        }

        let mut adb = db_paths
            .align_db
            .open(true)
            .context("failed to open align DB for merge")?;

        let mut adb_index = db_paths.align_db.with_extension("index").open(true)?;

        db_paths
            .align_db
            .with_extension("dbtype")
            .open(true)?
            .write_all(ALIGN_DBTYPE)?;

        let mut offset = 0;
        for prf_idx in 0..queries.len() {
            let mut n_written = 0;
            for prog_adb in &mut prog_adbs {
                let record_bytes = prog_adb.get(prf_idx)?;
                adb.write_all(record_bytes)?;
                n_written += record_bytes.len();
            }

            adb.write_all(&[0])?;
            n_written += 1;

            let s = format!("{}\t{offset}\t{}\n", prf_idx, n_written);
            adb_index.write_all(s.as_bytes())?;
            offset += n_written;
        }
    }

    // ---

    let align_tsv = match &args.io_args.seeds_output_path {
        Some(path) => path,
        None => &args.io_args.temp_dir_path.join("seeds.tsv"),
    };

    let now = Instant::now();
    run_mmseqs_convertalis(
        &db_paths.query_db,
        &db_paths.target_db,
        &db_paths.align_db,
        align_tsv,
        args,
    )
    .context("mmseqs convertalis failed")?;

    stats.set_mmseqs_time(crate::stats::MmseqsTimed::Convertalis, now.elapsed());

    // ---

    let now = Instant::now();
    let seeds = Seeds2::from_path(align_tsv).context("failed to build seeds")?;

    stats.set_mmseqs_time(crate::stats::MmseqsTimed::Index, now.elapsed());
    stats.set_mmseqs_time(crate::stats::MmseqsTimed::Total, time_start.elapsed());

    Ok(seeds)
}
