use crate::align::bounded::structs::bound::{join_bounds, CloudBoundGroup};
use crate::align::bounded::structs::row_bound_params::RowBoundParams;
use crate::align::bounded::structs::{CloudMatrixLinear, CloudSearchParams};
use crate::align::bounded::{
    backward_bounded, cloud_search_backward, cloud_search_forward, optimal_accuracy_bounded,
    posterior_bounded, traceback_bounded,
};
use crate::align::forward_bounded;
use crate::structs::dp_matrix::DpMatrix;
use crate::structs::{Alignment, DpMatrix3D, DpMatrixFlat, Profile, Sequence, Trace};
use crate::util::Average;
use anyhow::Result;
use std::fs::{create_dir_all, File};
use std::io::{stdout, BufWriter};
use std::time::Instant;

pub struct BoundedTimings {
    pub num_alignments: usize,
    pub reuse_times: Vec<usize>,
    pub cloud_forward_times: Vec<usize>,
    pub cloud_backward_times: Vec<usize>,
    pub forward_times: Vec<usize>,
    pub backward_times: Vec<usize>,
    pub posterior_times: Vec<usize>,
    pub optimal_times: Vec<usize>,
}

impl BoundedTimings {
    pub fn new(num_alignments: usize) -> Self {
        Self {
            num_alignments,
            reuse_times: vec![0; num_alignments],
            cloud_forward_times: vec![0; num_alignments],
            cloud_backward_times: vec![0; num_alignments],
            forward_times: vec![0; num_alignments],
            backward_times: vec![0; num_alignments],
            posterior_times: vec![0; num_alignments],
            optimal_times: vec![0; num_alignments],
        }
    }

    pub fn print_avg(&self) {
        println!("bounded pipeline timing summary");
        println!("-----------------------------");
        println!("reuse:          {:7}μs", self.reuse_times.avg());
        println!("cloud_forward:  {:7}μs", self.cloud_forward_times.avg());
        println!("cloud_backward: {:7}μs", self.cloud_backward_times.avg());
        println!("forward:        {:7}μs", self.forward_times.avg());
        println!("backward:       {:7}μs", self.backward_times.avg());
        println!("posterior:      {:7}μs", self.posterior_times.avg());
        println!("optimal:        {:7}μs", self.optimal_times.avg());
        println!(
            "all:       {:7}μs",
            self.reuse_times.avg()
                + self.cloud_forward_times.avg()
                + self.cloud_backward_times.avg()
                + self.forward_times.avg()
                + self.backward_times.avg()
                + self.posterior_times.avg()
                + self.optimal_times.avg()
        );
    }
}

pub fn pipeline_bounded(profiles: &mut [Profile], targets: &[Sequence]) -> Result<Vec<Alignment>> {
    let mut cloud_search_params = CloudSearchParams::default();
    let dump = true;

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
    let mut forward_matrix = DpMatrixFlat::new(max_target_length, max_profile_length);
    let mut backward_matrix = DpMatrixFlat::new(max_target_length, max_profile_length);
    let mut posterior_matrix = DpMatrixFlat::new(max_target_length, max_profile_length);
    let mut optimal_matrix = DpMatrixFlat::new(max_target_length, max_profile_length);

    // TODO: this needs to be implemented
    // let mut trace = Trace::default();

    let mut timings = BoundedTimings::new(profiles.len());
    let mut alignments: Vec<Alignment> = vec![];

    for (alignment_cnt, (profile, target)) in profiles.iter_mut().zip(targets.iter()).enumerate() {
        // for profile in profiles.iter_mut() {
        //     for target in targets.iter() {

        create_dir_all(format!("./out/{}", profile.name))?;

        profile.configure_for_target_length(target.length);

        let now = Instant::now();

        // TODO: this method might need work
        cloud_matrix.reuse(profile.length);

        forward_bounds.reuse(target.length, profile.length);
        backward_bounds.reuse(target.length, profile.length);

        forward_matrix.reuse(target.length, profile.length);
        backward_matrix.reuse(target.length, profile.length);
        posterior_matrix.reuse(target.length, profile.length);
        optimal_matrix.reuse(target.length, profile.length);

        timings.reuse_times[alignment_cnt] = now.elapsed().as_micros() as usize;

        // TODO: ***BIGTIME TODO***
        //       these need to be passed in on a vector to be set for each profile/target pair
        cloud_search_params.target_start = 1;
        cloud_search_params.profile_start = 1;
        cloud_search_params.target_end = target.length;
        cloud_search_params.profile_end = profile.length;

        let now = Instant::now();
        cloud_search_forward(
            profile,
            target,
            &mut cloud_matrix,
            &cloud_search_params,
            &mut forward_bounds,
        )?;
        timings.cloud_forward_times[alignment_cnt] = now.elapsed().as_micros() as usize;

        let now = Instant::now();
        cloud_search_backward(
            profile,
            target,
            &mut cloud_matrix,
            &cloud_search_params,
            &mut backward_bounds,
        )?;
        timings.cloud_backward_times[alignment_cnt] = now.elapsed().as_micros() as usize;

        join_bounds(&mut forward_bounds, &backward_bounds)?;

        forward_bounds.trim_wings();

        let row_bound_params = RowBoundParams::new(&forward_bounds);

        let now = Instant::now();
        forward_bounded(profile, target, &mut forward_matrix, &row_bound_params);
        timings.forward_times[alignment_cnt] = now.elapsed().as_micros() as usize;

        if dump {
            let mut forward_out = BufWriter::new(File::create(format!(
                "./out/{}/{}-fwd-bounded.mtx",
                profile.name, target.name
            ))?);
            forward_matrix.dump(&mut forward_out)?;
        }

        let now = Instant::now();
        backward_bounded(profile, target, &mut backward_matrix, &row_bound_params);
        timings.backward_times[alignment_cnt] = now.elapsed().as_micros() as usize;

        if dump {
            let mut backward_out = BufWriter::new(File::create(format!(
                "./out/{}/{}-bwd-bounded.mtx",
                profile.name, target.name
            ))?);
            backward_matrix.dump(&mut backward_out)?;
        }

        let now = Instant::now();
        posterior_bounded(
            profile,
            &forward_matrix,
            &backward_matrix,
            &mut posterior_matrix,
            &row_bound_params,
        );
        timings.posterior_times[alignment_cnt] = now.elapsed().as_micros() as usize;

        if dump {
            let mut posterior_out = BufWriter::new(File::create(format!(
                "./out/{}/{}-post-bounded.mtx",
                profile.name, target.name
            ))?);
            posterior_matrix.dump(&mut posterior_out)?;
        }

        let now = Instant::now();
        optimal_accuracy_bounded(
            profile,
            &posterior_matrix,
            &mut optimal_matrix,
            &row_bound_params,
        );
        timings.optimal_times[alignment_cnt] = now.elapsed().as_micros() as usize;

        if dump {
            let mut optimal_out = BufWriter::new(File::create(format!(
                "./out/{}/{}-opt-bounded.mtx",
                profile.name, target.name
            ))?);
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
            let mut trace_out = BufWriter::new(File::create(format!(
                "./out/{}/{}-trace-bounded.out",
                profile.name, target.name
            ))?);
            trace.dump(&mut trace_out, profile, target)?;
        }

        alignments.push(Alignment::new(&trace, profile, target));
    }
    timings.print_avg();
    Ok(alignments)
}
