use std::{fs::create_dir_all, io::Write, time::Instant};

use anyhow::Context;
use libnail::align::{structs::Seed, Bits};

use crate::{
    args::SearchArgs,
    io::{Database, Fasta, P7Hmm, Seeds2},
    mmseqs::{
        consts::{ALIGN_DBTYPE, PREFILTER_DBTYPE},
        run_mmseqs_align, run_mmseqs_convertalis, run_mmseqs_prefilter, run_mmseqs_search,
        write_mmseqs_profile_database, write_mmseqs_sequence_database, ByteBuffer, MmseqsDbPaths,
        MmseqsScoreModel, PrefilterDb,
    },
    pipeline::StageResult,
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

pub fn seed_profile_to_sequence_progressive(
    profiles: &P7Hmm,
    seqs: &Fasta,
    stats: &mut Stats,
    args: &SearchArgs,
) -> anyhow::Result<Seeds2> {
    let time_start = Instant::now();

    let paths = MmseqsDbPaths::new(&args.io_args.temp_dir_path);

    // ---

    write_mmseqs_sequence_database(seqs, &paths.target_db)?;
    write_mmseqs_profile_database(profiles.values(), &paths.query_db)?;

    let now = Instant::now();
    run_mmseqs_prefilter(
        &paths.query_db,
        &paths.target_db,
        &paths.prefilter_db,
        None,
        args,
    )?;
    stats.set_mmseqs_time(crate::stats::MmseqsTimed::Prefilter, now.elapsed());

    // ---

    let mut pf = PrefilterDb::from_path(&paths.prefilter_db)
        .context("failed to open initial prefilterDB")?;

    #[derive(PartialEq, Debug)]
    enum State {
        Active,
        Final(usize),
        Terminated,
    }

    let dir = &paths.prefilter_db.with_file_name("prog/");

    let mut i = 0;
    let mut n_take = 200;
    let mut state: Vec<State> = (0..profiles.len()).map(|_| State::Active).collect();
    let mut prog_adbs = vec![];
    while state.contains(&State::Active) {
        let pf_path = dir.join(format!("{i}/prefilterDB"));
        let dir = pf_path.parent().unwrap();
        create_dir_all(dir)?;

        {
            // note: scoped to drop file handles and force a write
            let mut prog_pfdb = pf_path.open(true)?;
            let mut prog_pfdb_index = pf_path.with_extension("index").open(true)?;
            pf_path
                .with_extension("dbtype")
                .open(true)?
                .write_all(PREFILTER_DBTYPE)?;

            let mut offset = 0;
            for (prf_idx, prf_state) in state.iter_mut().enumerate() {
                let record_bytes = match prf_state {
                    State::Active => match pf.next_n(prf_idx, n_take)? {
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

        let prog_adb_path = pf_path.with_file_name("alignDB");

        let now = Instant::now();
        run_mmseqs_align(
            &paths.query_db,
            &paths.target_db,
            &pf_path,
            &prog_adb_path,
            None,
            args,
        )?;
        stats.add_mmseqs_time(crate::stats::MmseqsTimed::Align, now.elapsed());

        let mut prog_adb =
            PrefilterDb::from_path(prog_adb_path).context("failed to open alignDB")?;
        let mut report_out = pf_path.with_file_name("report.txt").open(true)?;

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

                    if frac < 0.01 {
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
        let mut adb = paths.align_db.open(true)?;
        let mut adb_index = paths.align_db.with_extension("index").open(true)?;
        paths
            .align_db
            .with_extension("dbtype")
            .open(true)?
            .write_all(ALIGN_DBTYPE)?;

        let mut offset = 0;
        for prf_idx in 0..profiles.len() {
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
        None => &paths.dir()?.join("seeds.tsv"),
    };

    let now = Instant::now();
    run_mmseqs_convertalis(
        &paths.query_db,
        &paths.target_db,
        &paths.align_db,
        align_tsv,
        args,
    )?;
    stats.set_mmseqs_time(crate::stats::MmseqsTimed::Convertalis, now.elapsed());

    let now = Instant::now();
    let seeds = Seeds2::from_path(align_tsv).context("failed to build seeds")?;
    stats.set_mmseqs_time(crate::stats::MmseqsTimed::Index, now.elapsed());

    stats.set_mmseqs_time(crate::stats::MmseqsTimed::Total, time_start.elapsed());

    Ok(seeds)
}

pub fn seed_profile_to_sequence(
    profiles: &P7Hmm,
    seqs: &Fasta,
    args: &SearchArgs,
) -> anyhow::Result<Seeds2> {
    let paths = MmseqsDbPaths::new(&args.io_args.temp_dir_path);

    write_mmseqs_sequence_database(seqs, &paths.target_db)?;

    let align_tsv = match &args.io_args.seeds_output_path {
        Some(path) => path,
        None => &paths.dir()?.join("seeds.tsv"),
    };

    write_mmseqs_profile_database(profiles.values(), &paths.query_db)?;
    run_mmseqs_search(&paths, align_tsv, args, MmseqsScoreModel::Profile)?;
    let seeds = Seeds2::from_path(align_tsv)?;
    Ok(seeds)
}

pub fn seed_sequence_to_sequence(
    queries: &Fasta,
    targets: &Fasta,
    args: &SearchArgs,
) -> anyhow::Result<Seeds2> {
    let paths = MmseqsDbPaths::new(&args.io_args.temp_dir_path);

    write_mmseqs_sequence_database(targets, &paths.target_db)?;
    write_mmseqs_sequence_database(queries, &paths.query_db)?;

    let align_tsv = match &args.io_args.seeds_output_path {
        Some(path) => path,
        None => &paths.dir()?.join("seeds.tsv"),
    };

    run_mmseqs_search(&paths, align_tsv, args, MmseqsScoreModel::Blosum62)?;
    let seeds = Seeds2::from_path(align_tsv)?;
    Ok(seeds)
}
