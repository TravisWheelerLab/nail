use crate::extension_traits::PathBufExt;
use crate::pipeline::{AlignArgs, SeedMap};
use libnail::align::bounded::structs::{
    CloudBoundGroup, CloudMatrixLinear, CloudSearchParams, DpMatrixSparse, RowBounds, Seed,
};
use libnail::align::bounded::{
    backward_bounded, cloud_search_backward, cloud_search_forward, forward_bounded, null1_score,
    null2_score_bounded, optimal_accuracy_bounded, posterior_bounded, traceback_bounded,
};
use libnail::structs::alignment::ScoreParams;
use libnail::structs::{Alignment, Profile, Sequence, Trace};
use rayon::iter::{IntoParallelIterator, ParallelIterator};
use std::collections::HashMap;
use std::fs::File;
use std::io::{stdout, BufWriter, Write};
use std::sync::Mutex;

pub fn align_threaded_bounded(
    args: &AlignArgs,
    mut profiles: Vec<Profile>,
    targets: Vec<Sequence>,
    seed_map: SeedMap,
) -> anyhow::Result<()> {
    let tab_results_writer: Mutex<BufWriter<File>> =
        Mutex::new(args.output_args.tsv_results_path.open(true)?);
    let ali_results_writer: Mutex<Box<dyn Write + Send>> = match args.output_args.ali_results_path {
        Some(ref path) => Mutex::new(Box::new(path.open(true)?)),
        None => Mutex::new(Box::new(stdout())),
    };

    let dp = AlignmentStructs::default();

    let score_params = ScoreParams::new(targets.len());

    let cloud_search_params = CloudSearchParams {
        gamma: args.libnail_args.gamma,
        alpha: args.libnail_args.alpha,
        beta: args.libnail_args.beta,
    };

    let mut target_map: HashMap<String, Sequence> = HashMap::new();
    for target in targets {
        target_map.insert(target.name.clone(), target);
    }

    let mut profile_seeds_pairs: Vec<(&mut Profile, &Vec<Seed>)> = vec![];

    for profile in profiles.iter_mut() {
        match seed_map.get(&profile.name) {
            Some(seeds) => profile_seeds_pairs.push((profile, seeds)),
            None => {
                continue;
            }
        }
    }

    #[derive(Default, Clone)]
    struct AlignmentStructs {
        cloud_matrix: CloudMatrixLinear,
        forward_bounds: CloudBoundGroup,
        backward_bounds: CloudBoundGroup,
        forward_matrix: DpMatrixSparse,
        backward_matrix: DpMatrixSparse,
        posterior_matrix: DpMatrixSparse,
        optimal_matrix: DpMatrixSparse,
    }

    profile_seeds_pairs.into_par_iter().for_each_with(
        (dp, score_params),
        |(dp, score_params), (profile, seeds)| {
            for seed in seeds {
                let target = target_map.get(&seed.target_name).unwrap();
                profile.configure_for_target_length(target.length);

                dp.cloud_matrix.reuse(profile.length);
                dp.forward_bounds.reuse(target.length, profile.length);
                dp.backward_bounds.reuse(target.length, profile.length);

                cloud_search_forward(
                    profile,
                    target,
                    seed,
                    &mut dp.cloud_matrix,
                    &cloud_search_params,
                    &mut dp.forward_bounds,
                );

                cloud_search_backward(
                    profile,
                    target,
                    seed,
                    &mut dp.cloud_matrix,
                    &cloud_search_params,
                    &mut dp.backward_bounds,
                );

                CloudBoundGroup::join_bounds(&mut dp.forward_bounds, &dp.backward_bounds);

                if !dp.forward_bounds.valid() {
                    println!("cloud bound fail");
                    continue;
                }

                dp.forward_bounds.trim_wings();

                let row_bounds = RowBounds::new(&dp.forward_bounds);

                if !row_bounds.valid() {
                    println!("row bound fail");
                    continue;
                }

                dp.forward_matrix
                    .reuse(target.length, profile.length, &row_bounds);
                dp.backward_matrix
                    .reuse(target.length, profile.length, &row_bounds);
                dp.posterior_matrix
                    .reuse(target.length, profile.length, &row_bounds);
                dp.optimal_matrix
                    .reuse(target.length, profile.length, &row_bounds);

                // we use the forward score to compute the final bit score (later)
                score_params.forward_score_nats =
                    forward_bounded(profile, target, &mut dp.forward_matrix, &row_bounds);

                backward_bounded(profile, target, &mut dp.backward_matrix, &row_bounds);

                posterior_bounded(
                    profile,
                    &dp.forward_matrix,
                    &dp.backward_matrix,
                    &mut dp.posterior_matrix,
                    &row_bounds,
                );

                score_params.null_score_nats = null1_score(target.length);
                score_params.bias_correction_score_nats =
                    null2_score_bounded(&dp.posterior_matrix, profile, target, &row_bounds);

                optimal_accuracy_bounded(
                    profile,
                    &dp.posterior_matrix,
                    &mut dp.optimal_matrix,
                    &row_bounds,
                );

                let mut trace = Trace::new(target.length, profile.length);
                traceback_bounded(
                    profile,
                    &dp.posterior_matrix,
                    &dp.optimal_matrix,
                    &mut trace,
                    row_bounds.target_end,
                );

                let alignment = Alignment::from_trace(&trace, profile, target, score_params);

                if alignment.evalue <= args.output_args.evalue_threshold {
                    match tab_results_writer.lock() {
                        Ok(mut writer) => writeln!(writer, "{}", alignment.tab_string())
                            .expect("failed to write tabular output"),
                        Err(_) => panic!("tabular results writer mutex was poisoned"),
                    };

                    match ali_results_writer.lock() {
                        Ok(mut writer) => {
                            writeln!(writer, "{}", alignment.ali_string())
                                .expect("failed to write alignment output");
                        }
                        Err(_) => panic!("alignment results writer mutex was poisoned"),
                    }
                }
            }
        },
    );

    Ok(())
}
