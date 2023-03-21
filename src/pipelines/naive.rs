use crate::align::bounded::traceback_bounded;
use crate::align::{backward, forward, optimal_accuracy, posterior, traceback};
use crate::structs::dp_matrix::DpMatrix;
use crate::structs::{Alignment, DpMatrix3D, DpMatrixFlat, Profile, Sequence, Trace};
use crate::util::Average;
use anyhow::Result;
use std::fs::{create_dir_all, File};
use std::io::{stdout, BufWriter};
use std::time::Instant;

pub struct NaiveTimings {
    pub num_alignments: usize,
    pub reuse_times: Vec<usize>,
    pub forward_times: Vec<usize>,
    pub backward_times: Vec<usize>,
    pub posterior_times: Vec<usize>,
    pub optimal_times: Vec<usize>,
}

impl NaiveTimings {
    pub fn new(num_alignments: usize) -> Self {
        Self {
            num_alignments,
            reuse_times: vec![0; num_alignments],
            forward_times: vec![0; num_alignments],
            backward_times: vec![0; num_alignments],
            posterior_times: vec![0; num_alignments],
            optimal_times: vec![0; num_alignments],
        }
    }

    pub fn print_avg(&self) {
        println!("naive pipeline timing summary");
        println!("-----------------------------");
        println!("reuse:     {:7}μs", self.reuse_times.avg());
        println!("forward:   {:7}μs", self.forward_times.avg());
        println!("backward:  {:7}μs", self.backward_times.avg());
        println!("posterior: {:7}μs", self.posterior_times.avg());
        println!("optimal:   {:7}μs", self.optimal_times.avg());
        println!(
            "all:       {:7}μs",
            self.reuse_times.avg()
                + self.forward_times.avg()
                + self.backward_times.avg()
                + self.posterior_times.avg()
                + self.optimal_times.avg()
        );
    }
}

pub fn pipeline_naive(profiles: &mut [Profile], targets: &[Sequence]) -> Result<Vec<Alignment>> {
    let dump = true;

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

    let mut timings = NaiveTimings::new(profiles.len());
    let mut alignments: Vec<Alignment> = vec![];

    for (alignment_cnt, (profile, target)) in profiles.iter_mut().zip(targets.iter()).enumerate() {
        // for profile in profiles.iter_mut() {
        //     for target in targets.iter() {

        create_dir_all(format!("./out/{}", profile.name))?;

        profile.configure_for_target_length(target.length);

        let now = Instant::now();

        forward_matrix.reuse(target.length, profile.length);
        backward_matrix.reuse(target.length, profile.length);
        posterior_matrix.reuse(target.length, profile.length);
        optimal_matrix.reuse(target.length, profile.length);

        timings.reuse_times[alignment_cnt] = now.elapsed().as_micros() as usize;

        let now = Instant::now();
        forward(profile, target, &mut forward_matrix);
        timings.forward_times[alignment_cnt] = now.elapsed().as_micros() as usize;

        if dump {
            let mut forward_out = BufWriter::new(File::create(format!(
                "./out/{}/{}-fwd-naive.mtx",
                profile.name, target.name
            ))?);
            forward_matrix.dump(&mut forward_out)?;
        }

        let now = Instant::now();
        backward(profile, target, &mut backward_matrix);
        timings.backward_times[alignment_cnt] = now.elapsed().as_micros() as usize;

        if dump {
            let mut backward_out = BufWriter::new(File::create(format!(
                "./out/{}/{}-bwd-naive.mtx",
                profile.name, target.name
            ))?);
            backward_matrix.dump(&mut backward_out)?;
        }

        let now = Instant::now();
        posterior(
            profile,
            &forward_matrix,
            &backward_matrix,
            &mut posterior_matrix,
        );
        timings.posterior_times[alignment_cnt] = now.elapsed().as_micros() as usize;

        if dump {
            let mut posterior_out = BufWriter::new(File::create(format!(
                "./out/{}/{}-post-naive.mtx",
                profile.name, target.name
            ))?);
            posterior_matrix.dump(&mut posterior_out)?;
        }

        let now = Instant::now();
        optimal_accuracy(profile, &posterior_matrix, &mut optimal_matrix);
        timings.optimal_times[alignment_cnt] = now.elapsed().as_micros() as usize;

        if dump {
            let mut optimal_out = BufWriter::new(File::create(format!(
                "./out/{}/{}-opt-naive.mtx",
                profile.name, target.name
            ))?);
            optimal_matrix.dump(&mut optimal_out)?;
        }

        let mut trace = Trace::new(target.length, profile.length);
        // traceback(profile, &posterior_matrix, &optimal_matrix, &mut trace);
        traceback_bounded(
            profile,
            &posterior_matrix,
            &optimal_matrix,
            &mut trace,
            target.length,
        );

        if dump {
            let mut trace_out = BufWriter::new(File::create(format!(
                "./out/{}/{}-trace-naive.out",
                profile.name, target.name
            ))?);
            trace.dump(&mut trace_out, profile, target)?;
        }

        alignments.push(Alignment::new(&trace, profile, target));
    }

    timings.print_avg();
    Ok(alignments)
}
