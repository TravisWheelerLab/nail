use crate::align::{backward, forward, optimal_accuracy, posterior, traceback};
use crate::structs::{Alignment, DpMatrix, Profile, Sequence, Trace};
use anyhow::Result;
use std::fs::File;
use std::io::{stdout, BufWriter};
use std::time::Instant;

pub fn pipeline_naive(profiles: &mut [Profile], targets: &[Sequence]) -> Result<()> {
    let dump = false;

    let max_profile_length = profiles
        .iter()
        .fold(0usize, |acc: usize, p: &Profile| acc.max(p.length));

    let max_target_length = targets
        .iter()
        .fold(0usize, |acc: usize, s: &Sequence| acc.max(s.length));

    let mut forward_matrix = DpMatrix::new(max_target_length, max_profile_length);
    let mut backward_matrix = DpMatrix::new(max_target_length, max_profile_length);
    let mut posterior_matrix = DpMatrix::new(max_target_length, max_profile_length);
    let mut optimal_matrix = DpMatrix::new(max_target_length, max_profile_length);

    for (profile, target) in profiles.iter_mut().zip(targets.iter()) {
        // for profile in profiles.iter_mut() {
        //     for target in targets.iter() {
        profile.configure_for_target_length(target.length);

        forward_matrix.reuse(target.length, profile.length);
        backward_matrix.reuse(target.length, profile.length);
        posterior_matrix.reuse(target.length, profile.length);
        optimal_matrix.reuse(target.length, profile.length);

        forward(profile, target, &mut forward_matrix);

        if dump {
            let mut forward_out = BufWriter::new(File::create("./nale-dump/forward.mtx")?);
            forward_matrix.dump(&mut forward_out)?;
        }

        backward(profile, target, &mut backward_matrix);

        if dump {
            let mut backward_out = BufWriter::new(File::create("./nale-dump/backward.mtx")?);
            backward_matrix.dump(&mut backward_out)?;
        }

        posterior(
            profile,
            &forward_matrix,
            &backward_matrix,
            &mut posterior_matrix,
        );

        if dump {
            let mut posterior_out = BufWriter::new(File::create("./nale-dump/posterior.mtx")?);
            posterior_matrix.dump(&mut posterior_out)?;
        }

        optimal_accuracy(profile, &posterior_matrix, &mut optimal_matrix);

        if dump {
            let mut optimal_out = BufWriter::new(File::create("./nale-dump/optimal.mtx")?);
            optimal_matrix.dump(&mut optimal_out)?;
        }

        let mut trace = Trace::new(target.length, profile.length);
        traceback(profile, &posterior_matrix, &optimal_matrix, &mut trace);

        if dump {
            let mut trace_out = BufWriter::new(File::create("./nale-dump/trace.dump")?);
            trace.dump(&mut trace_out, profile, target)?;
        }

        let alignment = Alignment::new(&trace, profile, target);
        // alignment.dump(&mut stdout())?;
        // }
    }

    Ok(())
}
