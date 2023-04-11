use crate::align::bounded::structs::{
    CloudBoundGroup, CloudDebugAnnotations, CloudMatrixLinear, CloudSearchParams, RowBounds,
};
use crate::align::bounded::{
    backward_bounded, cloud_search_backward, cloud_search_forward, forward_bounded,
    optimal_accuracy_bounded, posterior_bounded, traceback_bounded,
};
use crate::output::path_buf_ext::PathBufExt;
use crate::structs::{Alignment, DpMatrixFlat, Profile, Sequence, Trace};
use anyhow::Result;
use std::collections::HashMap;
use std::fmt::{Display, Formatter};
use std::io::Write;
use std::path::PathBuf;

pub struct BoundedPipelineParams {
    pub write_debug: bool,
    pub allow_overwrite: bool,
    pub debug_path: PathBuf,
    pub cloud_search_params: CloudSearchParams,
}

impl Default for BoundedPipelineParams {
    fn default() -> Self {
        Self {
            write_debug: false,
            allow_overwrite: false,
            debug_path: PathBuf::from("./nale-debug"),
            cloud_search_params: CloudSearchParams::default(),
        }
    }
}

// TODO: rename and put this elsewhere
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

    let mut cloud_debug = if params.write_debug {
        Some(CloudDebugAnnotations::default())
    } else {
        None
    };

    for profile_name in &profile_names {
        let profile = profile_map.get_mut(*profile_name).unwrap();
        let seeds = profile_seeds.get(*profile_name).unwrap();
        for seed in seeds {
            let target = target_map.get(&seed.target_name[..]).unwrap();

            println!("{}", profile.name);
            println!("{}", seed);

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

            if let Some(ref mut cloud_debug) = cloud_debug {
                cloud_debug.add_forward_bounds(&forward_bounds);
                cloud_debug.add_backward_bounds(&backward_bounds);
            }

            CloudBoundGroup::join_bounds(&mut forward_bounds, &backward_bounds)?;

            forward_bounds.trim_wings();

            let row_bounds = RowBounds::new(&forward_bounds);

            if let Some(ref mut cloud_debug) = cloud_debug {
                cloud_debug.add_joined_bounds(&forward_bounds);
                cloud_debug.add_row_bounds(&row_bounds);

                let mut out = params
                    .debug_path
                    .join(&format!("{}-{}.bounds.json", profile_name, target.name))
                    .open(params.allow_overwrite)?;

                let serialized = serde_json::to_string(&cloud_debug).unwrap();
                write!(out, "{serialized}")?;
            }

            forward_matrix.reuse(target.length, profile.length);
            backward_matrix.reuse(target.length, profile.length);
            posterior_matrix.reuse(target.length, profile.length);
            optimal_matrix.reuse(target.length, profile.length);

            // forward_matrix.reuse(target.length, profile.length, &row_bound_params);
            // backward_matrix.reuse(target.length, profile.length, &row_bound_params);
            // posterior_matrix.reuse(target.length, profile.length, &row_bound_params);
            // optimal_matrix.reuse(target.length, profile.length, &row_bound_params);

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

            alignments.push(Alignment::new(&trace, profile, target));
        }
    }
    Ok(alignments)
}
