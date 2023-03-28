use crate::align::bounded::traceback_bounded;
use crate::align::naive::backward::backward;
use crate::align::naive::forward::forward;
use crate::align::naive::optimal_accuracy::optimal_accuracy;
use crate::align::naive::posterior::posterior;
use crate::output::output_debug::{
    get_profile_target_output_dir_path, set_file_name_and_get_buf_writer,
};
use crate::structs::dp_matrix::DpMatrix;
use crate::structs::{Alignment, DpMatrixFlat, Profile, Sequence, Trace};
use anyhow::Result;
use std::path::PathBuf;

pub struct NaivePipelineParams {
    pub write_debug: bool,
    pub allow_overwrite: bool,
    pub root_debug_dir_path: PathBuf,
}

impl Default for NaivePipelineParams {
    fn default() -> Self {
        Self {
            write_debug: false,
            allow_overwrite: false,
            root_debug_dir_path: PathBuf::from("./nale-debug"),
        }
    }
}

pub fn pipeline_naive(
    profiles: &mut [Profile],
    targets: &[Sequence],
    params: &NaivePipelineParams,
) -> Result<Vec<Alignment>> {
    let max_profile_length = profiles
        .iter()
        .fold(0usize, |acc: usize, p: &Profile| acc.max(p.length));

    let max_target_length = targets
        .iter()
        .fold(0usize, |acc: usize, s: &Sequence| acc.max(s.length));

    let mut forward_matrix = DpMatrixFlat::new(max_target_length, max_profile_length);
    let mut backward_matrix = DpMatrixFlat::new(max_target_length, max_profile_length);
    let mut posterior_matrix = DpMatrixFlat::new(max_target_length, max_profile_length);
    let mut optimal_matrix = DpMatrixFlat::new(max_target_length, max_profile_length);

    let mut alignments: Vec<Alignment> = vec![];

    for (profile, target) in profiles.iter_mut().zip(targets.iter()) {
        // for profile in profiles.iter_mut() {
        //     for target in targets.iter() {

        let mut current_debug_dir_path = if params.write_debug {
            let path = get_profile_target_output_dir_path(
                &params.root_debug_dir_path,
                &profile.name,
                &target.name,
            )?;
            Some(path)
        } else {
            None
        };

        profile.configure_for_target_length(target.length);

        forward_matrix.reuse(target.length, profile.length);
        backward_matrix.reuse(target.length, profile.length);
        posterior_matrix.reuse(target.length, profile.length);
        optimal_matrix.reuse(target.length, profile.length);

        forward(profile, target, &mut forward_matrix)?;

        backward(profile, target, &mut backward_matrix)?;

        if params.write_debug {
            if let Some(ref mut path) = current_debug_dir_path {
                let out = &mut set_file_name_and_get_buf_writer(
                    path,
                    "full-forward.mtx",
                    params.allow_overwrite,
                )?;
                forward_matrix.dump(out)?;

                let out = &mut set_file_name_and_get_buf_writer(
                    path,
                    "full-backward.mtx",
                    params.allow_overwrite,
                )?;
                backward_matrix.dump(out)?;
            }
        }

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

        if params.write_debug {
            if let Some(ref mut path) = current_debug_dir_path {
                let out = &mut set_file_name_and_get_buf_writer(
                    path,
                    "full-posterior.mtx",
                    params.allow_overwrite,
                )?;
                posterior_matrix.dump(out)?;

                let out = &mut set_file_name_and_get_buf_writer(
                    path,
                    "full-optimal.mtx",
                    params.allow_overwrite,
                )?;
                optimal_matrix.dump(out)?;

                let out = &mut set_file_name_and_get_buf_writer(
                    path,
                    "full-trace.mtx",
                    params.allow_overwrite,
                )?;
                trace.dump(out, profile, target)?;
            }
        }

        alignments.push(Alignment::new(&trace, profile, target));
    }

    Ok(alignments)
}
