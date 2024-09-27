use std::time::Instant;

use libnail::{
    align::{
        backward, forward, null_one_score, null_two_score, optimal_accuracy, p_value, posterior,
        structs::{Alignment, AlignmentBuilder, DpMatrixSparse, RowBounds, Trace},
        traceback,
    },
    structs::{Profile, Sequence},
};

use crate::args::SearchArgs;

#[derive(Clone)]
pub struct AlignConfig {
    pub do_null_two: bool,
}

impl Default for AlignConfig {
    fn default() -> Self {
        Self { do_null_two: true }
    }
}

pub trait AlignStep: dyn_clone::DynClone {
    fn run(
        &mut self,
        profile: &mut Profile,
        target: &Sequence,
        bounds: &RowBounds,
    ) -> anyhow::Result<Alignment>;
}

dyn_clone::clone_trait_object!(AlignStep);

#[derive(Default, Clone)]
pub struct DefaultAlignStep {
    forward_matrix: DpMatrixSparse,
    backward_matrix: DpMatrixSparse,
    posterior_matrix: DpMatrixSparse,
    optimal_matrix: DpMatrixSparse,
    forward_p_value_threshold: f64,
    target_count: usize,
    e_value_threshold: f64,
    config: AlignConfig,
}

impl DefaultAlignStep {
    pub fn new(args: &SearchArgs) -> Self {
        Self {
            target_count: match args.nail_args.target_database_size {
                Some(size) => size,
                None => panic!(),
            },
            forward_p_value_threshold: args.nail_args.forward_pvalue_threshold,
            e_value_threshold: args.output_args.evalue_threshold,
            ..Default::default()
        }
    }
}

impl AlignStep for DefaultAlignStep {
    fn run(
        &mut self,
        profile: &mut Profile,
        target: &Sequence,
        bounds: &RowBounds,
    ) -> anyhow::Result<Alignment> {
        // configuring for the target length
        // adjusts special state transitions
        profile.configure_for_target_length(target.length);

        let now = Instant::now();
        self.forward_matrix
            .reuse(target.length, profile.length, bounds);
        let mut init_time = now.elapsed();

        // we use the forward score to compute the final bit score (later)
        let now = Instant::now();
        let forward_score = forward(profile, target, &mut self.forward_matrix, bounds).to_bits()
            // the denominator is the null one score
            - null_one_score(target.length);
        let forward_time = now.elapsed();

        // for now we compute the P-value for filtering purposes
        let forward_pvalue = p_value(forward_score, profile.forward_lambda, profile.forward_tau);

        let alignment = AlignmentBuilder::default()
            .with_profile(profile)
            .with_target(target)
            .with_database_size(self.target_count)
            .with_cell_count(bounds.num_cells)
            .with_forward_score(forward_score)
            .with_init_time(init_time)
            .with_forward_time(forward_time);

        if forward_pvalue >= self.forward_p_value_threshold {
            return alignment.build();
        }

        let now = Instant::now();
        self.backward_matrix
            .reuse(target.length, profile.length, bounds);
        self.posterior_matrix
            .reuse(target.length, profile.length, bounds);
        self.optimal_matrix
            .reuse(target.length, profile.length, bounds);
        init_time += now.elapsed();

        let now = Instant::now();
        backward(profile, target, &mut self.backward_matrix, bounds);
        let backward_time = now.elapsed();

        let now = Instant::now();
        posterior(
            profile,
            &self.forward_matrix,
            &self.backward_matrix,
            &mut self.posterior_matrix,
            bounds,
        );
        let posterior_time = now.elapsed();

        let now = Instant::now();
        optimal_accuracy(
            profile,
            &self.posterior_matrix,
            &mut self.optimal_matrix,
            bounds,
        );
        let optimal_accuracy_time = now.elapsed();

        let now = Instant::now();
        let mut trace = Trace::new(target.length, profile.length);
        traceback(
            profile,
            &self.posterior_matrix,
            &self.optimal_matrix,
            &mut trace,
            bounds.target_end,
        );
        let traceback_time = now.elapsed();

        let now = Instant::now();
        let null_two = if self.config.do_null_two {
            Some(null_two_score(
                &self.posterior_matrix,
                profile,
                target,
                bounds,
            ))
        } else {
            None
        };
        let null_two_time = now.elapsed();

        alignment
            .with_trace(&trace)
            .with_null_two(null_two)
            .with_init_time(init_time)
            .with_backward_time(backward_time)
            .with_posterior_time(posterior_time)
            .with_optimal_accuracy_time(optimal_accuracy_time)
            .with_traceback_time(traceback_time)
            .with_null_two_time(null_two_time)
            .build()
    }
}
