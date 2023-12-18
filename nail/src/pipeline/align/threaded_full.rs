use crate::extension_traits::PathBufExt;
use crate::pipeline::{AlignArgs, SeedMap};
use libnail::align::bounded::structs::Seed;
use libnail::align::bounded::{null1_score, traceback_bounded};
use libnail::align::naive::backward::backward;
use libnail::align::naive::forward::forward;
use libnail::align::naive::null2::null2_score;
use libnail::align::naive::optimal_accuracy::optimal_accuracy;
use libnail::align::naive::posterior::posterior;
use libnail::structs::alignment::ScoreParams;
use libnail::structs::{Alignment, DpMatrixFlat, Profile, Sequence, Trace};
use rayon::iter::{IntoParallelIterator, ParallelIterator};
use std::collections::HashMap;
use std::fs::File;
use std::io::{stdout, BufWriter, Write};
use std::sync::Mutex;

pub fn align_threaded_full(
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
        forward_matrix: DpMatrixFlat,
        backward_matrix: DpMatrixFlat,
        posterior_matrix: DpMatrixFlat,
        optimal_matrix: DpMatrixFlat,
    }

    profile_seeds_pairs.into_par_iter().for_each_with(
        (dp, score_params),
        |(dp, score_params), (profile, seeds)| {
            for seed in seeds {
                let target = target_map.get(&seed.target_name).unwrap();
                profile.configure_for_target_length(target.length);

                dp.forward_matrix.reuse(target.length, profile.length);
                dp.backward_matrix.reuse(target.length, profile.length);
                dp.posterior_matrix.reuse(target.length, profile.length);
                dp.optimal_matrix.reuse(target.length, profile.length);

                // we use the forward score to compute the final bit score (later)
                score_params.forward_score_nats = forward(profile, target, &mut dp.forward_matrix);

                backward(profile, target, &mut dp.backward_matrix);

                posterior(
                    profile,
                    &dp.forward_matrix,
                    &dp.backward_matrix,
                    &mut dp.posterior_matrix,
                );

                score_params.null_score_nats = null1_score(target.length);
                score_params.bias_correction_score_nats =
                    null2_score(&dp.posterior_matrix, profile, target);

                optimal_accuracy(profile, &dp.posterior_matrix, &mut dp.optimal_matrix);

                let mut trace = Trace::new(target.length, profile.length);
                traceback_bounded(
                    profile,
                    &dp.posterior_matrix,
                    &dp.optimal_matrix,
                    &mut trace,
                    target.length,
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
