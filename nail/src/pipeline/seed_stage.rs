use std::{
    fs::{create_dir_all, File},
    io::{BufWriter, Write},
};

use anyhow::Context;
use libnail::structs::Profile;

use crate::{
    args::SearchArgs,
    io::{Database, Fasta, P7Hmm, SeedList, Seeds},
    mmseqs::{
        consts::PREFILTER_DBTYPE, run_mmseqs_align, run_mmseqs_prefilter, run_mmseqs_search,
        write_mmseqs_profile_database, write_mmseqs_sequence_database, ByteBuffer, MmseqsDbPaths,
        MmseqsScoreModel, PrefilterDb,
    },
    util::PathBufExt,
};

pub fn seed_profile_to_sequence_progressive(
    profiles: &P7Hmm,
    seqs: &Fasta,
    args: &SearchArgs,
) -> anyhow::Result<Seeds> {
    let paths = MmseqsDbPaths::new(&args.io_args.temp_dir_path);

    // ---

    write_mmseqs_sequence_database(seqs, &paths.target_db)?;
    write_mmseqs_profile_database(profiles.values(), &paths.query_db)?;

    run_mmseqs_prefilter(
        &paths.query_db,
        &paths.target_db,
        &paths.prefilter_db,
        None,
        args,
    )?;

    // ---

    let mut pf = PrefilterDb::from_path(&paths.prefilter_db)
        .context("failed to open initial prefilterDB")?;

    #[derive(PartialEq, Debug)]
    enum State {
        Active,
        Final(usize),
        Terminated,
    }

    let mut state: Vec<State> = (0..pf.len()).map(|_| State::Active).collect();

    let dir = &paths.prefilter_db.with_file_name("prog/");

    let mut i = 0;
    let mut n_take = 200;
    while state.contains(&State::Active) {
        let start = std::time::Instant::now();

        let pf_path = dir.join(format!("{i}/prefilterDB"));
        let pf_index_path = pf_path.with_extension("index");
        let pf_dbtype_path = pf_path.with_extension("dbtype");

        let dir = pf_path.parent().unwrap();
        create_dir_all(dir)?;

        {
            // note: scoped to drop file handles and force a write
            let mut db_out = pf_path.open(true)?;
            let mut idx_out = pf_index_path.open(true)?;
            pf_dbtype_path.open(true)?.write_all(PREFILTER_DBTYPE)?;

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

                db_out.write_all(record_bytes)?;
                // seed groups are terminated with a null byte
                db_out.write_all(&[0])?;

                let n_written = record_bytes.len() + 1;

                let s = format!("{}\t{offset}\t{}\n", prf_idx, n_written);
                idx_out.write_all(s.as_bytes())?;

                offset += n_written;
            }
        }

        let align_tsv_path = pf_path.with_file_name("align.tsv");
        let align_db_path = pf_path.with_file_name("alignDB");

        let now = std::time::Instant::now();
        run_mmseqs_align(
            &paths.query_db,
            &paths.target_db,
            &pf_path,
            &align_db_path,
            align_tsv_path,
            None,
            args,
        )?;
        let align_time = now.elapsed();

        let mut adb = PrefilterDb::from_path(align_db_path).context("failed to open alignDB")?;
        let mut out = pf_path.with_file_name("report.txt").open(true)?;

        for (prf_idx, prf_state) in state.iter_mut().enumerate() {
            match prf_state {
                State::Active => {
                    let record_bytes = adb.get(prf_idx)?;
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
                    writeln!(out, "{prf_idx}: {cnt} / {n_take} | {frac:0.3}")?;
                }
                State::Final(n) => {
                    writeln!(out, "{prf_idx}: ({n})")?;
                    *prf_state = State::Terminated;
                }
                State::Terminated => continue,
            }
        }

        writeln!(out, "total: {:?}", start.elapsed())?;
        writeln!(out, "align: {:?}", align_time)?;

        i += 1;
        n_take *= 2;
    }

    // ---

    let align_tsv = match &args.io_args.seeds_output_path {
        Some(path) => path,
        None => &paths.dir()?.join("seeds.tsv"),
    };

    let tmp_tsv = align_tsv.with_file_name("seeds.tmp.tsv");

    let mut seeds_writer = BufWriter::new(File::create(align_tsv)?);

    let seeds = Seeds::from_path(align_tsv)?;

    Ok(seeds)
}

pub fn seed_profile_to_sequence(
    profiles: &P7Hmm,
    seqs: &Fasta,
    args: &SearchArgs,
) -> anyhow::Result<Seeds> {
    let paths = MmseqsDbPaths::new(&args.io_args.temp_dir_path);

    write_mmseqs_sequence_database(seqs, &paths.target_db)?;

    let align_tsv = match &args.io_args.seeds_output_path {
        Some(path) => path,
        None => &paths.dir()?.join("seeds.tsv"),
    };

    write_mmseqs_profile_database(profiles.values(), &paths.query_db)?;
    run_mmseqs_search(&paths, align_tsv, args, MmseqsScoreModel::Profile)?;
    let seeds = Seeds::from_path(align_tsv)?;
    Ok(seeds)
}

pub fn seed_sequence_to_sequence(
    queries: &Fasta,
    targets: &Fasta,
    args: &SearchArgs,
) -> anyhow::Result<Seeds> {
    let paths = MmseqsDbPaths::new(&args.io_args.temp_dir_path);

    write_mmseqs_sequence_database(targets, &paths.target_db)?;
    write_mmseqs_sequence_database(queries, &paths.query_db)?;

    let align_tsv = match &args.io_args.seeds_output_path {
        Some(path) => path,
        None => &paths.dir()?.join("seeds.tsv"),
    };

    run_mmseqs_search(&paths, align_tsv, args, MmseqsScoreModel::Blosum62)?;
    let seeds = Seeds::from_path(align_tsv)?;
    Ok(seeds)
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
