use crate::align::{backward, forward, optimal_accuracy, posterior, traceback};
use crate::structs::{Alignment, DpMatrix, Profile, Sequence, Trace};
use anyhow::Result;
use std::fs::File;
use std::io::{stdout, BufWriter};
use std::time::Instant;

pub fn pipeline_naive(profiles: &mut [Profile], targets: &[Sequence]) -> Result<()> {
    for profile in profiles.iter_mut() {
        for target in targets.iter() {
            profile.configure_for_length(target.length);

            // TODO: reuse matrices

            let mut forward_matrix = DpMatrix::new(target.length, profile.length);
            forward(profile, target, &mut forward_matrix);

            // let mut forward_out = BufWriter::new(File::create("./nale-dump/forward.mtx")?);
            // forward_matrix.dump(&mut forward_out)?;

            let mut backward_matrix = DpMatrix::new(target.length, profile.length);
            backward(profile, target, &mut backward_matrix);

            // let mut backward_out = BufWriter::new(File::create("./nale-dump/backward.mtx")?);
            // backward_matrix.dump(&mut backward_out)?;

            let mut posterior_matrix = DpMatrix::new(target.length, profile.length);
            posterior(
                profile,
                &forward_matrix,
                &backward_matrix,
                &mut posterior_matrix,
            );

            // let mut posterior_out = BufWriter::new(File::create("./nale-dump/posterior.mtx")?);
            // posterior_matrix.dump(&mut posterior_out)?;

            let mut optimal_matrix = DpMatrix::new(target.length, profile.length);
            optimal_accuracy(profile, &posterior_matrix, &mut optimal_matrix);

            // let mut optimal_out = BufWriter::new(File::create("./nale-dump/optimal.mtx")?);
            // optimal_matrix.dump(&mut optimal_out)?;

            let mut trace = Trace::new(profile.length, target.length);
            traceback(profile, &posterior_matrix, &optimal_matrix, &mut trace);

            let mut trace_out = BufWriter::new(File::create("./nale-dump/trace.dump")?);
            trace.dump(&mut trace_out, profile, target)?;

            // let alignment = Alignment::new(&trace, profile, target);
            // alignment.dump(&mut stdout())?;
        }
    }

    Ok(())
}
