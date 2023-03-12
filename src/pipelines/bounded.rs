use crate::align::bounded::structs::row_bound_params::RowBoundParams;
use crate::align::bounded::structs::{CloudMatrixLinear, CloudSearchParams};
use crate::align::bounded::{
    backward_bounded, cloud_search, optimal_accuracy_bounded, posterior_bounded,
};
use crate::align::{forward_bounded, traceback};
use crate::structs::{Alignment, DpMatrix, Profile, Sequence, Trace};
use anyhow::Result;
use std::fs::File;
use std::io::{stdout, BufWriter};
use std::time::Instant;

pub fn pipeline_bounded(profiles: &mut [Profile], targets: &[Sequence]) -> Result<()> {
    let mut cloud_search_params = CloudSearchParams::default();

    for profile in profiles.iter_mut() {
        for target in targets.iter() {
            profile.configure_for_length(target.length);

            cloud_search_params.target_start = 1;
            cloud_search_params.profile_start = 1;
            cloud_search_params.target_end = target.length;
            cloud_search_params.profile_end = profile.length;

            let mut cloud_matrix = CloudMatrixLinear::new(profile.length);

            let joined_bounds =
                cloud_search(profile, target, &mut cloud_matrix, &cloud_search_params)?;

            let row_bound_params = RowBoundParams::new(&joined_bounds);

            let mut forward_matrix = DpMatrix::new(profile.length, target.length);
            forward_bounded(profile, target, &mut forward_matrix, &row_bound_params);

            // let mut forward_out = BufWriter::new(File::create("./nale-bounded-dump/forward.mtx")?);
            // forward_matrix.dump(&mut forward_out)?;

            let mut backward_matrix = DpMatrix::new(profile.length, target.length);
            backward_bounded(profile, target, &mut backward_matrix, &row_bound_params);

            // let mut backward_out =
            //     BufWriter::new(File::create("./nale-bounded-dump/backward.mtx")?);
            // backward_matrix.dump(&mut backward_out)?;

            let mut posterior_matrix = DpMatrix::new(target.length, profile.length);
            posterior_bounded(
                profile,
                &forward_matrix,
                &backward_matrix,
                &mut posterior_matrix,
                &row_bound_params,
            );

            // let mut posterior_out =
            //     BufWriter::new(File::create("./nale-bounded-dump/posterior.mtx")?);
            // posterior_matrix.dump(&mut posterior_out)?;

            let mut optimal_matrix = DpMatrix::new(target.length, profile.length);
            optimal_accuracy_bounded(
                profile,
                &posterior_matrix,
                &mut optimal_matrix,
                &row_bound_params,
            );

            // let mut optimal_out = BufWriter::new(File::create("./nale-bounded-dump/optimal.mtx")?);
            // optimal_matrix.dump(&mut optimal_out)?;

            let mut trace = Trace::new(profile.length, target.length);
            traceback(profile, &posterior_matrix, &optimal_matrix, &mut trace);

            // let mut trace_out = BufWriter::new(File::create("./nale-bounded-dump/trace.dump")?);
            // trace.dump(&mut trace_out, profile, target)?;

            // let alignment = Alignment::new(&trace, profile, target);
            // alignment.dump(&mut stdout())?;
        }
    }
    Ok(())
}
