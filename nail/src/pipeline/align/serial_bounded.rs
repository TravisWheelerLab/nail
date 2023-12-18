use crate::extension_traits::PathBufExt;
use crate::pipeline::{AlignArgs, SeedMap, TargetNotFoundError};
use libnail::align::bounded::structs::{
    CloudBoundGroup, CloudMatrixLinear, CloudSearchParams, DpMatrixSparse, RowBounds,
};
use libnail::align::bounded::{
    backward_bounded, cloud_search_backward, cloud_search_forward, forward_bounded, null1_score,
    null2_score_bounded, optimal_accuracy_bounded, posterior_bounded, traceback_bounded,
};
use libnail::structs::alignment::ScoreParams;
use libnail::structs::{Alignment, Profile, Sequence, Trace};
use std::collections::HashMap;
use std::io::{stdout, Write};

pub fn align_serial_bounded(
    args: &AlignArgs,
    mut profiles: Vec<Profile>,
    targets: Vec<Sequence>,
    seed_map: SeedMap,
) -> anyhow::Result<()> {
    let cloud_search_params = CloudSearchParams {
        gamma: args.libnail_args.gamma,
        alpha: args.libnail_args.alpha,
        beta: args.libnail_args.beta,
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

    let mut cloud_matrix = CloudMatrixLinear::new(max_profile_length);

    let mut forward_bounds = CloudBoundGroup::new(max_target_length, max_profile_length);
    let mut backward_bounds = CloudBoundGroup::new(max_target_length, max_profile_length);

    let mut forward_matrix =
        DpMatrixSparse::new(max_target_length, max_profile_length, &RowBounds::default());
    let mut backward_matrix =
        DpMatrixSparse::new(max_target_length, max_profile_length, &RowBounds::default());
    let mut posterior_matrix =
        DpMatrixSparse::new(max_target_length, max_profile_length, &RowBounds::default());
    let mut optimal_matrix =
        DpMatrixSparse::new(max_target_length, max_profile_length, &RowBounds::default());

    let mut tab_results_writer = args.output_args.tsv_results_path.open(true)?;

    let mut ali_results_writer: Box<dyn Write> = match args.output_args.ali_results_path {
        Some(ref path) => Box::new(path.open(true)?),
        None => Box::new(stdout()),
    };

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

            cloud_matrix.reuse(profile.length);
            forward_bounds.reuse(target.length, profile.length);
            backward_bounds.reuse(target.length, profile.length);

            cloud_search_forward(
                profile,
                target,
                seed,
                &mut cloud_matrix,
                &cloud_search_params,
                &mut forward_bounds,
            );

            cloud_search_backward(
                profile,
                target,
                seed,
                &mut cloud_matrix,
                &cloud_search_params,
                &mut backward_bounds,
            );

            CloudBoundGroup::join_bounds(&mut forward_bounds, &backward_bounds);

            if !forward_bounds.valid() {
                println!("cloud bound fail");
                continue;
            }

            forward_bounds.trim_wings();

            let row_bounds = RowBounds::new(&forward_bounds);

            if !row_bounds.valid() {
                println!("row bound fail");
                continue;
            }

            forward_matrix.reuse(target.length, profile.length, &row_bounds);
            backward_matrix.reuse(target.length, profile.length, &row_bounds);
            posterior_matrix.reuse(target.length, profile.length, &row_bounds);
            optimal_matrix.reuse(target.length, profile.length, &row_bounds);

            // we use the forward score to compute the final bit score (later)
            score_params.forward_score_nats =
                forward_bounded(profile, target, &mut forward_matrix, &row_bounds);

            backward_bounded(profile, target, &mut backward_matrix, &row_bounds);

            posterior_bounded(
                profile,
                &forward_matrix,
                &backward_matrix,
                &mut posterior_matrix,
                &row_bounds,
            );

            optimal_accuracy_bounded(profile, &posterior_matrix, &mut optimal_matrix, &row_bounds);

            let mut trace = Trace::new(target.length, profile.length);
            traceback_bounded(
                profile,
                &posterior_matrix,
                &optimal_matrix,
                &mut trace,
                row_bounds.target_end,
            );

            score_params.null_score_nats = null1_score(target.length);
            score_params.bias_correction_score_nats =
                null2_score_bounded(&posterior_matrix, profile, target, &row_bounds);

            let alignment = Alignment::from_trace(&trace, profile, target, &score_params);

            if alignment.evalue <= args.output_args.evalue_threshold {
                writeln!(ali_results_writer, "{}", alignment.ali_string())?;
                writeln!(tab_results_writer, "{}", alignment.tab_string())?;
            }
        }
    }
    Ok(())
}
