use std::{cell::RefCell, collections::HashMap};

use anyhow::Context;
use libnail::{
    align::{
        backward, cloud_search_bwd, cloud_search_fwd, forward, posterior,
        structs::{AdMatrixLinear, Cloud, DpMatrix, DpMatrixSparse, Relationship, RowBounds},
        CloudSearchParams,
    },
    structs::Profile,
};
use rayon::{iter::ParallelIterator, slice::ParallelSlice};
use regex::Regex;
use thread_local::ThreadLocal;

use crate::{
    args::SearchArgs,
    io::Fasta,
    pipeline::Pipeline,
    search::{build_pipeline, read_queries, seed, Queries},
    stats::Stats,
    util::PathExt,
};

pub fn dev_search(mut args: SearchArgs) -> anyhow::Result<()> {
    let queries = read_queries(&args.query_path)?;
    let targets = Fasta::from_path(&args.target_path).context("failed to read target fasta")?;

    let mut stats = Stats::new(queries.len(), targets.len());

    match args.expert_args.target_database_size {
        Some(_) => {}
        None => args.expert_args.target_database_size = Some(targets.len()),
    }

    let seeds = seed(&queries, &targets, &mut stats, &mut args)?;

    if args.pipeline_args.only_seed {
        return Ok(());
    }

    let pipeline = build_pipeline(queries, targets, stats, &mut args)?;
    let tl_pipeline: ThreadLocal<RefCell<Pipeline>> = ThreadLocal::new();

    seeds
        .seeds
        .par_chunks(100)
        .panic_fuse()
        .try_for_each(|chunk| -> anyhow::Result<()> {
            let mut pipeline = tl_pipeline
                .get_or(|| RefCell::new(pipeline.clone()))
                .borrow_mut();

            let res = chunk
                .iter()
                .map(|seed| pipeline.run(seed))
                .collect::<Result<Vec<_>, _>>()?;

            pipeline.output.run(&res)?;

            Ok(())
        })?;

    Ok(())
}

pub fn dev_mx(mut args: SearchArgs) -> anyhow::Result<()> {
    let queries = read_queries(&args.query_path)?;
    let mut targets = Fasta::from_path(&args.target_path).context("failed to read target fasta")?;

    let mut stats = Stats::new(queries.len(), targets.len());

    match args.expert_args.target_database_size {
        Some(_) => {}
        None => args.expert_args.target_database_size = Some(targets.len()),
    }

    let seeds = seed(&queries, &targets, &mut stats, &mut args)?;

    let mut profiles: HashMap<String, Profile> = match queries {
        Queries::Sequence(fasta) => fasta
            .par_iter()
            .filter_map(|s| Profile::from_blosum_62_and_seq(&s).ok())
            .collect::<Vec<_>>(),
        Queries::Profile(p7hmm) => p7hmm.par_iter().collect::<Vec<_>>(),
    }
    .into_iter()
    .map(|p| (p.name.clone(), p))
    .collect();

    let mx_dir = std::path::PathBuf::from("./mx");
    mx_dir.create_dir()?;

    let mut cld_mx = AdMatrixLinear::default();
    let mut fwd_cld = Cloud::default();
    let mut bwd_cld = Cloud::default();

    let params = CloudSearchParams {
        gamma: args.pipeline_args.gamma,
        alpha: args.pipeline_args.alpha,
        beta: args.pipeline_args.beta,
    };

    let mut fwd_mx = DpMatrixSparse::default();
    let mut bwd_mx = DpMatrixSparse::default();
    let mut post_mx = DpMatrixSparse::default();
    let mut opt_mx = DpMatrixSparse::default();

    let file_re = Regex::new(r"[^A-Za-z0-9._-]").unwrap();
    '_main: for seed in seeds.seeds.iter() {
        let prf = profiles.get_mut(&seed.prf).unwrap();
        let seq = targets.get(&seed.seq).unwrap();

        let prefix = mx_dir.join(format!(
            "{}-{}",
            file_re.replace_all(&prf.name, "_"),
            file_re.replace_all(&seq.name, "_")
        ));

        let mut row_bounds = RowBounds::new(seq.length);

        if args.dev_args.full_dp {
            row_bounds.fill_rectangle(1, 1, seq.length, prf.length);
        } else {
            cld_mx.reuse(seq.length);
            fwd_cld.reuse(seq.length, prf.length);
            bwd_cld.reuse(seq.length, prf.length);

            for attempts in 0..5 {
                cld_mx.reuse(seq.length);
                fwd_cld.reuse(seq.length, prf.length);
                bwd_cld.reuse(seq.length, prf.length);

                let params = params.scale(1.0 + (attempts as f32 * 0.5));
                let _ = cloud_search_fwd(prf, &seq, seed, &mut cld_mx, &params, &mut fwd_cld);

                cld_mx.reuse(seq.length);

                let _ = cloud_search_bwd(prf, &seq, seed, &mut cld_mx, &params, &mut bwd_cld);

                match fwd_cld.anti_diagonal_relationship(&bwd_cld) {
                    Relationship::Disjoint(_) => {
                        if attempts >= 4 {
                            continue '_main;
                        }
                    }
                    Relationship::Intersecting(_) => break,
                };
            }

            fwd_cld.merge(&bwd_cld);
            fwd_cld.square_corners();
            let trim_result = fwd_cld.trim_wings();

            match trim_result {
                Ok(_) => {
                    row_bounds.fill_from_cloud(&fwd_cld);
                }
                Err(_) => {
                    panic!("row bounds failure");
                }
            }
        }

        fwd_mx.reuse(seq.length, prf.length, &row_bounds);
        bwd_mx.reuse(seq.length, prf.length, &row_bounds);
        post_mx.reuse(seq.length, prf.length, &row_bounds);
        opt_mx.reuse(seq.length, prf.length, &row_bounds);

        let _ = forward(prf, &seq, &mut fwd_mx, &row_bounds);
        backward(prf, &seq, &mut bwd_mx, &row_bounds);
        posterior(prf, &fwd_mx, &bwd_mx, &mut post_mx, &row_bounds);

        fwd_mx.dump(&mut prefix.with_extension("fwd.mx").open(true)?)?;
        bwd_mx.dump(&mut prefix.with_extension("bwd.mx").open(true)?)?;
        post_mx.dump(&mut prefix.with_extension("post.mx").open(true)?)?;
    }

    Ok(())
}
