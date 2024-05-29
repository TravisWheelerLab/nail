use libnail::{
    align::{
        backward, forward, null_one_score, null_two_score, optimal_accuracy, p_value, posterior,
        structs::{Alignment, AlignmentBuilder, DpMatrixSparse, RowBounds, Trace},
        traceback,
    },
    structs::{Profile, Sequence},
};

use crate::args::AlignArgs;

pub trait AlignStep: dyn_clone::DynClone {
    fn run(
        &mut self,
        profile: &mut Profile,
        target: &Sequence,
        bounds: &RowBounds,
    ) -> Option<Alignment>;
}

dyn_clone::clone_trait_object!(AlignStep);

#[derive(Default, Clone)]
pub struct DefaultAlignStep {
    forward_matrix: DpMatrixSparse,
    backward_matrix: DpMatrixSparse,
    posterior_matrix: DpMatrixSparse,
    optimal_matrix: DpMatrixSparse,
    forward_pvalue_threshold: f64,
    target_count: usize,
    e_value_threshold: f64,
}

impl DefaultAlignStep {
    pub fn new(args: &AlignArgs, target_count: usize) -> Self {
        Self {
            target_count: match args.nail_args.target_database_size {
                Some(size) => size,
                None => target_count,
            },
            forward_pvalue_threshold: args.nail_args.forward_pvalue_threshold,
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
    ) -> Option<Alignment> {
        // configuring for the target length
        // adjusts special state transitions
        profile.configure_for_target_length(target.length);

        self.forward_matrix
            .reuse(target.length, profile.length, bounds);
        self.backward_matrix
            .reuse(target.length, profile.length, bounds);
        self.posterior_matrix
            .reuse(target.length, profile.length, bounds);
        self.optimal_matrix
            .reuse(target.length, profile.length, bounds);

        // we use the forward score to compute the final bit score (later)
        let forward_score = forward(profile, target, &mut self.forward_matrix, bounds).to_bits();

        // for now we compute the P-value for filtering purposes
        let forward_pvalue = p_value(forward_score, profile.forward_lambda, profile.forward_tau);

        if forward_pvalue >= self.forward_pvalue_threshold {
            return None;
        }

        backward(profile, target, &mut self.backward_matrix, bounds);

        posterior(
            profile,
            &self.forward_matrix,
            &self.backward_matrix,
            &mut self.posterior_matrix,
            bounds,
        );

        optimal_accuracy(
            profile,
            &self.posterior_matrix,
            &mut self.optimal_matrix,
            bounds,
        );

        let mut trace = Trace::new(target.length, profile.length);
        traceback(
            profile,
            &self.posterior_matrix,
            &self.optimal_matrix,
            &mut trace,
            bounds.target_end,
        );

        let null_one = null_one_score(target.length);
        let null_two = null_two_score(&self.posterior_matrix, profile, target, bounds);

        let cell_fraction = bounds.num_cells() as f32 / (profile.length * target.length) as f32;

        let alignment = AlignmentBuilder::new(&trace)
            .with_profile(profile)
            .with_target(target)
            .with_target_count(self.target_count)
            .with_forward_score(forward_score)
            .with_null_one(null_one)
            .with_null_two(null_two)
            .with_cell_fraction(cell_fraction)
            .build()
            .unwrap();

        match alignment.e_value {
            Some(e_value) if e_value <= self.e_value_threshold => Some(alignment),
            _ => None,
        }
    }
}
