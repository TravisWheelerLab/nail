use crate::align::bounded::structs::{
    CloudBoundGroup, CloudMatrixLinear, CloudSearchParams, RowBoundParams,
};
use crate::align::bounded::{
    backward_bounded, cloud_search_backward, cloud_search_forward, forward_bounded,
    optimal_accuracy_bounded, posterior_bounded, traceback_bounded,
};
use crate::output::output_debug::{
    get_profile_target_output_dir_path, set_file_name_and_get_buf_writer,
};
use crate::structs::dp_matrix::DpMatrix;
use crate::structs::{Alignment, DpMatrixFlat, Profile, Sequence, Trace};
use crate::viz::SodaJson;
use anyhow::Result;
use std::collections::HashMap;
use std::fmt::{Display, Formatter};
use std::path::PathBuf;

pub struct BoundedPipelineParams {
    pub write_debug: bool,
    pub allow_overwrite: bool,
    pub root_debug_dir_path: PathBuf,
    pub cloud_search_params: CloudSearchParams,
}

impl Default for BoundedPipelineParams {
    fn default() -> Self {
        Self {
            write_debug: false,
            allow_overwrite: false,
            root_debug_dir_path: PathBuf::from("./nale-debug"),
            cloud_search_params: CloudSearchParams::default(),
        }
    }
}

pub struct Seed {
    pub target_name: String,
    pub target_start: usize,
    pub target_end: usize,
    pub profile_start: usize,
    pub profile_end: usize,
}

impl Display for Seed {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        writeln!(f, "{}", self.target_name)?;
        writeln!(f, "  {}, {}", self.target_start, self.profile_start)?;
        writeln!(f, "  {}, {}", self.target_end, self.profile_end)?;
        Ok(())
    }
}

pub fn pipeline_bounded(
    profile_map: &mut HashMap<String, Profile>,
    target_map: &HashMap<String, Sequence>,
    profile_seeds: &HashMap<String, Vec<Seed>>,
    params: &BoundedPipelineParams,
) -> Result<Vec<Alignment>> {
    let max_profile_length = profile_map
        .values()
        .fold(0usize, |acc: usize, p: &Profile| acc.max(p.length));

    let max_target_length = target_map
        .values()
        .fold(0usize, |acc: usize, s: &Sequence| acc.max(s.length));

    let mut cloud_matrix = CloudMatrixLinear::new(max_profile_length);

    let mut forward_bounds = CloudBoundGroup::new(max_target_length, max_profile_length);
    let mut backward_bounds = CloudBoundGroup::new(max_target_length, max_profile_length);

    // TODO: we probably want to overwrite these matrices
    //       i.e. forward becomes posterior, backward becomes optimal or something
    let mut forward_matrix = DpMatrixFlat::new(max_target_length, max_profile_length);
    let mut backward_matrix = DpMatrixFlat::new(max_target_length, max_profile_length);
    let mut posterior_matrix = DpMatrixFlat::new(max_target_length, max_profile_length);
    let mut optimal_matrix = DpMatrixFlat::new(max_target_length, max_profile_length);

    // let mut forward_matrix = DpMatrixSparse::default();
    // let mut backward_matrix = DpMatrixSparse::default();
    // let mut posterior_matrix = DpMatrixSparse::default();
    // let mut optimal_matrix = DpMatrixSparse::default();

    // TODO: this needs to be implemented
    // let mut trace = Trace::default();

    let mut alignments: Vec<Alignment> = vec![];

    let mut profile_names: Vec<&String> = profile_seeds.keys().collect();
    profile_names.sort();

    // for (i, profile_name) in profile_seeds.keys().enumerate() {
    for profile_name in &profile_names {
        let profile = profile_map.get_mut(*profile_name).unwrap();
        let seeds = profile_seeds.get(*profile_name).unwrap();
        for seed in seeds {
            let target = target_map.get(&seed.target_name[..]).unwrap();

            println!("{}", profile.name);
            println!("{}", seed);

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

            // TODO: this method might need work
            cloud_matrix.reuse(profile.length);

            forward_bounds.reuse(target.length, profile.length);
            backward_bounds.reuse(target.length, profile.length);

            cloud_search_forward(
                profile,
                target,
                seed,
                &mut cloud_matrix,
                &params.cloud_search_params,
                &mut forward_bounds,
            )?;

            cloud_search_backward(
                profile,
                target,
                seed,
                &mut cloud_matrix,
                &params.cloud_search_params,
                &mut backward_bounds,
            )?;

            if params.write_debug {
                if let Some(ref mut path) = current_debug_dir_path {
                    let out = &mut set_file_name_and_get_buf_writer(
                        path,
                        "forward-bounds.json",
                        params.allow_overwrite,
                    )?;
                    forward_bounds.soda_json(out)?;

                    let out = &mut set_file_name_and_get_buf_writer(
                        path,
                        "backward-bounds.json",
                        params.allow_overwrite,
                    )?;
                    backward_bounds.soda_json(out)?;
                }
            }

            CloudBoundGroup::join_bounds(&mut forward_bounds, &backward_bounds)?;

            forward_bounds.trim_wings();

            let row_bound_params = RowBoundParams::new(&forward_bounds);

            if params.write_debug {
                if let Some(ref mut path) = current_debug_dir_path {
                    let out = &mut set_file_name_and_get_buf_writer(
                        path,
                        "joined-bounds.json",
                        params.allow_overwrite,
                    )?;
                    forward_bounds.soda_json(out)?;

                    let out = &mut set_file_name_and_get_buf_writer(
                        path,
                        "row-bounds.txt",
                        params.allow_overwrite,
                    )?;
                    row_bound_params.dump(out)?;
                }
            }

            // if !row_bound_params.valid() {
            //     continue;
            // }

            forward_matrix.reuse(target.length, profile.length);
            backward_matrix.reuse(target.length, profile.length);
            posterior_matrix.reuse(target.length, profile.length);
            optimal_matrix.reuse(target.length, profile.length);

            // forward_matrix.reuse(target.length, profile.length, &row_bound_params);
            // backward_matrix.reuse(target.length, profile.length, &row_bound_params);
            // posterior_matrix.reuse(target.length, profile.length, &row_bound_params);
            // optimal_matrix.reuse(target.length, profile.length, &row_bound_params);

            forward_bounded(profile, target, &mut forward_matrix, &row_bound_params);

            backward_bounded(profile, target, &mut backward_matrix, &row_bound_params);

            if params.write_debug {
                if let Some(ref mut path) = current_debug_dir_path {
                    let out = &mut set_file_name_and_get_buf_writer(
                        path,
                        "bounded-forward.mtx",
                        params.allow_overwrite,
                    )?;
                    forward_matrix.dump(out)?;

                    let out = &mut set_file_name_and_get_buf_writer(
                        path,
                        "bounded-backward.mtx",
                        params.allow_overwrite,
                    )?;
                    backward_matrix.dump(out)?;
                }
            }

            posterior_bounded(
                profile,
                &forward_matrix,
                &backward_matrix,
                &mut posterior_matrix,
                &row_bound_params,
            );

            optimal_accuracy_bounded(
                profile,
                &posterior_matrix,
                &mut optimal_matrix,
                &row_bound_params,
            );

            let mut trace = Trace::new(target.length, profile.length);
            traceback_bounded(
                profile,
                &posterior_matrix,
                &optimal_matrix,
                &mut trace,
                row_bound_params.target_end,
            );

            if params.write_debug {
                if let Some(ref mut path) = current_debug_dir_path {
                    let out = &mut set_file_name_and_get_buf_writer(
                        path,
                        "bounded-posterior.mtx",
                        params.allow_overwrite,
                    )?;
                    posterior_matrix.dump(out)?;

                    let out = &mut set_file_name_and_get_buf_writer(
                        path,
                        "bounded-optimal.mtx",
                        params.allow_overwrite,
                    )?;
                    optimal_matrix.dump(out)?;

                    let out = &mut set_file_name_and_get_buf_writer(
                        path,
                        "bounded-trace.mtx",
                        params.allow_overwrite,
                    )?;
                    trace.dump(out, profile, target)?;
                }
            }

            alignments.push(Alignment::new(&trace, profile, target));
        }
    }
    Ok(alignments)
}
