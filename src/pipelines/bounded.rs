use crate::align::bounded::structs::bound::{join_bounds, CloudBoundGroup};
use crate::align::bounded::structs::row_bound_params::RowBoundParams;
use crate::align::bounded::structs::{CloudMatrixLinear, CloudSearchParams};
use crate::align::bounded::{
    backward_bounded, cloud_search_backward, cloud_search_forward, optimal_accuracy_bounded,
    posterior_bounded, traceback_bounded,
};
use crate::align::forward_bounded;
use crate::structs::{Alignment, DpMatrix, Profile, Sequence, Trace};
use anyhow::Result;
use std::fs::File;
use std::io::{stdout, BufWriter};
use std::time::Instant;

pub fn pipeline_bounded(profiles: &mut [Profile], targets: &[Sequence]) -> Result<()> {
    let mut cloud_search_params = CloudSearchParams::default();
    let dump = false;

    let max_profile_length = profiles
        .iter()
        .fold(0usize, |acc: usize, p: &Profile| acc.max(p.length));

    let max_target_length = targets
        .iter()
        .fold(0usize, |acc: usize, s: &Sequence| acc.max(s.length));

    let mut cloud_matrix = CloudMatrixLinear::new(max_profile_length);

    let mut forward_bounds = CloudBoundGroup::new(max_target_length, max_profile_length);
    let mut backward_bounds = CloudBoundGroup::new(max_target_length, max_profile_length);

    // TODO: (todo after sparse DpMatrix is implemented)
    //       we probably want to overwrite these matrices
    //       i.e. forward becomes posterior, backward becomes optimal or something
    let mut forward_matrix = DpMatrix::new(max_target_length, max_profile_length);
    let mut backward_matrix = DpMatrix::new(max_target_length, max_profile_length);
    let mut posterior_matrix = DpMatrix::new(max_target_length, max_profile_length);
    let mut optimal_matrix = DpMatrix::new(max_target_length, max_profile_length);
    // TODO: this needs to be implemented
    // let mut trace = Trace::default();

    for (profile, target) in profiles.iter_mut().zip(targets.iter()) {
        // println!("t: {}, p: {}", target.length, profile.length);

        // for profile in profiles.iter_mut() {
        //     for target in targets.iter() {
        profile.configure_for_target_length(target.length);

        // TODO: this method might need work
        cloud_matrix.reuse(profile.length);

        forward_bounds.reuse(target.length, profile.length);
        backward_bounds.reuse(target.length, profile.length);

        forward_matrix.reuse(target.length, profile.length);
        backward_matrix.reuse(target.length, profile.length);
        posterior_matrix.reuse(target.length, profile.length);
        optimal_matrix.reuse(target.length, profile.length);

        // TODO: ***BIGTIME TODO***
        //       these need to be passed in on a vector to be set for each profile/target pair
        cloud_search_params.target_start = 1;
        cloud_search_params.profile_start = 1;
        cloud_search_params.target_end = target.length;
        cloud_search_params.profile_end = profile.length;

        cloud_search_forward(
            profile,
            target,
            &mut cloud_matrix,
            &cloud_search_params,
            &mut forward_bounds,
        )?;

        cloud_search_backward(
            profile,
            target,
            &mut cloud_matrix,
            &cloud_search_params,
            &mut backward_bounds,
        )?;

        join_bounds(&mut forward_bounds, &backward_bounds)?;

        forward_bounds.trim_wings();

        let row_bound_params = RowBoundParams::new(&forward_bounds);

        forward_bounded(profile, target, &mut forward_matrix, &row_bound_params);

        if dump {
            let mut forward_out = BufWriter::new(File::create("./nale-bounded-dump/forward.mtx")?);
            forward_matrix.dump(&mut forward_out)?;
        }

        backward_bounded(profile, target, &mut backward_matrix, &row_bound_params);

        if dump {
            let mut backward_out =
                BufWriter::new(File::create("./nale-bounded-dump/backward.mtx")?);
            backward_matrix.dump(&mut backward_out)?;
        }

        posterior_bounded(
            profile,
            &forward_matrix,
            &backward_matrix,
            &mut posterior_matrix,
            &row_bound_params,
        );

        if dump {
            let mut posterior_out =
                BufWriter::new(File::create("./nale-bounded-dump/posterior.mtx")?);
            posterior_matrix.dump(&mut posterior_out)?;
        }

        optimal_accuracy_bounded(
            profile,
            &posterior_matrix,
            &mut optimal_matrix,
            &row_bound_params,
        );

        if dump {
            let mut optimal_out = BufWriter::new(File::create("./nale-bounded-dump/optimal.mtx")?);
            optimal_matrix.dump(&mut optimal_out)?;
        }

        let mut trace = Trace::new(target.length, profile.length);
        traceback_bounded(
            profile,
            &posterior_matrix,
            &optimal_matrix,
            &mut trace,
            row_bound_params.target_end,
        );

        if dump {
            let mut trace_out = BufWriter::new(File::create("./nale-bounded-dump/trace.dump")?);
            trace.dump(&mut trace_out, profile, target)?;
        }

        let alignment = Alignment::new(&trace, profile, target);
        // alignment.dump(&mut stdout())?;
        // }
    }

    Ok(())
}
