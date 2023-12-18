use crate::extension_traits::PathBufExt;
use crate::pipeline::{AlignArgs, SeedMap, TargetNotFoundError};
use libnail::align::bounded::{null1_score, traceback_bounded};
use libnail::align::naive::backward::backward;
use libnail::align::naive::forward::forward;
use libnail::align::naive::null2::null2_score;
use libnail::align::naive::optimal_accuracy::optimal_accuracy;
use libnail::align::naive::posterior::posterior;
use libnail::structs::alignment::ScoreParams;
use libnail::structs::{Alignment, DpMatrixFlat, Profile, Sequence, Trace};
use std::collections::HashMap;
use std::io::{stdout, Write};

pub fn align_serial_full(
    args: &AlignArgs,
    mut profiles: Vec<Profile>,
    targets: Vec<Sequence>,
    seed_map: SeedMap,
) -> anyhow::Result<()> {
    let mut tab_results_writer = args.output_args.tsv_results_path.open(true)?;

    let mut ali_results_writer: Box<dyn Write> = match args.output_args.ali_results_path {
        Some(ref path) => Box::new(path.open(true)?),
        None => Box::new(stdout()),
    };

    let mut score_params = ScoreParams::new(targets.len());

    let mut target_map: HashMap<String, Sequence> = HashMap::new();
    for target in targets {
        target_map.insert(target.name.clone(), target);
    }

    let max_profile_length = profiles
        .iter()
        .fold(0usize, |acc: usize, p: &Profile| acc.max(p.length));

    let max_target_length = target_map
        .values()
        .fold(0usize, |acc: usize, s: &Sequence| acc.max(s.length));

    let mut forward_matrix = DpMatrixFlat::new(max_target_length, max_profile_length);
    let mut backward_matrix = DpMatrixFlat::new(max_target_length, max_profile_length);
    let mut posterior_matrix = DpMatrixFlat::new(max_target_length, max_profile_length);
    let mut optimal_matrix = DpMatrixFlat::new(max_target_length, max_profile_length);

    for profile in profiles.iter_mut() {
        let seeds = match seed_map.get(&profile.name) {
            Some(seeds) => seeds,
            None => {
                continue;
            }
        };

        for seed in seeds {
            let target =
                target_map
                    .get(&seed.target_name[..])
                    .ok_or_else(|| TargetNotFoundError {
                        target_name: seed.target_name.clone(),
                    })?;

            profile.configure_for_target_length(target.length);

            forward_matrix.reuse(target.length, profile.length);
            backward_matrix.reuse(target.length, profile.length);
            posterior_matrix.reuse(target.length, profile.length);
            optimal_matrix.reuse(target.length, profile.length);

            // we use the forward score to compute the final bit score (later)
            score_params.forward_score_nats = forward(profile, target, &mut forward_matrix);

            backward(profile, target, &mut backward_matrix);

            posterior(
                profile,
                &forward_matrix,
                &backward_matrix,
                &mut posterior_matrix,
            );

            optimal_accuracy(profile, &posterior_matrix, &mut optimal_matrix);

            let mut trace = Trace::new(target.length, profile.length);
            traceback_bounded(
                profile,
                &posterior_matrix,
                &optimal_matrix,
                &mut trace,
                target.length,
            );

            score_params.null_score_nats = null1_score(target.length);
            score_params.bias_correction_score_nats =
                null2_score(&posterior_matrix, profile, target);

            let alignment = Alignment::from_trace(&trace, profile, target, &score_params);

            if alignment.evalue <= args.output_args.evalue_threshold {
                writeln!(ali_results_writer, "{}", alignment.ali_string())?;
                writeln!(tab_results_writer, "{}", alignment.tab_string())?;
            }
        }
    }
    Ok(())
}
