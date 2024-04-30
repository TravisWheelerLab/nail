mod align;
use std::collections::HashMap;

pub use align::*;
mod prep;
use libnail::{
    align::{
        backward, cloud_search_backward, cloud_search_forward, composition_bias_score, forward,
        length_bias_score, optimal_accuracy, posterior,
        structs::{AntiDiagonalBounds, CloudMatrixLinear, DpMatrixSparse, RowBounds, Seed, Trace},
        traceback, CloudSearchParams, CloudSearchScores,
    },
    structs::{Profile, Sequence},
};
pub use prep::*;
mod search;
pub use search::*;
mod seed;
pub use seed::*;

struct SeedStep {
    seeds: HashMap<String, HashMap<String, Seed>>,
}

impl SeedStep {
    fn run(&self, profile: &Profile, target: &Sequence) -> Option<Seed> {
        Some(self.seeds.get(&profile.name)?.get(&target.name)?.clone())
    }
}

struct CloudSearchStep {
    cloud_matrix: CloudMatrixLinear,
    forward_bounds: AntiDiagonalBounds,
    reverse_bounds: AntiDiagonalBounds,
    row_bounds: RowBounds,
    params: CloudSearchParams,
}

impl CloudSearchStep {
    // TODO: move this to libnail
    fn cloud_score(forward_scores: &CloudSearchScores, reverse_scores: &CloudSearchScores) -> f32 {
        // this approximates the score for the forward
        // cloud that extends past the seed end point
        let disjoint_forward_score = forward_scores.max_score - forward_scores.max_score_within;

        // this approximates the score for the reverse
        // cloud that extends past the seed start point
        let disjoint_reverse_score = reverse_scores.max_score - reverse_scores.max_score_within;

        // this approximates the score of the intersection
        // of the forward and reverse clouds
        let intersection_score = forward_scores
            .max_score_within
            .max(reverse_scores.max_score_within);

        let cloud_score_nats = intersection_score + disjoint_forward_score + disjoint_reverse_score;

        // convert to bits
        cloud_score_nats / std::f32::consts::LN_2
    }

    fn run(&mut self, profile: &Profile, target: &Sequence, seed: &Seed) -> Option<&RowBounds> {
        self.cloud_matrix.reuse(profile.length);
        self.forward_bounds.reuse(target.length, profile.length);
        self.reverse_bounds.reuse(target.length, profile.length);
        self.row_bounds.reuse(target.length);

        let forward_scores = cloud_search_forward(
            profile,
            target,
            seed,
            &mut self.cloud_matrix,
            &self.params,
            &mut self.forward_bounds,
        );

        let reverse_scores = cloud_search_backward(
            profile,
            target,
            seed,
            &mut self.cloud_matrix,
            &self.params,
            &mut self.reverse_bounds,
        );

        // TODO: make this actually work
        // if Self::cloud_score(&forward_scores, &reverse_scores) < self.params.threshold {
        //     return None;
        // }

        let bounds_intersect =
            self.forward_bounds.max_anti_diagonal_idx >= self.reverse_bounds.min_anti_diagonal_idx;

        if bounds_intersect {
            AntiDiagonalBounds::join_merge(&mut self.forward_bounds, &self.reverse_bounds);
            if !self.forward_bounds.valid() {
                self.forward_bounds.fill_rectangle(
                    seed.target_start,
                    seed.profile_start,
                    seed.target_end,
                    seed.profile_end,
                );
            }

            self.forward_bounds.square_corners();
            self.forward_bounds.trim_wings();

            self.row_bounds
                .fill_from_anti_diagonal_bounds(&self.forward_bounds);

            if !self.row_bounds.valid() {
                self.row_bounds.fill_rectangle(
                    seed.target_start,
                    seed.profile_start,
                    seed.target_end,
                    seed.profile_end,
                );
            }
        } else {
            self.row_bounds.fill_rectangle(
                seed.target_start,
                seed.profile_start,
                seed.target_end,
                seed.profile_end,
            );
        }

        Some(&self.row_bounds)
    }
}

struct AlignmentStep {
    forward_matrix: DpMatrixSparse,
    backward_matrix: DpMatrixSparse,
    posterior_matrix: DpMatrixSparse,
    optimal_matrix: DpMatrixSparse,
}

impl AlignmentStep {
    fn run(&mut self, profile: &Profile, target: &Sequence, bounds: &RowBounds) {
        todo!();

        // TODO: need to refactor profile to produce a transition score struct
        // profile.configure_for_target_length(target.length);

        // TODO: move filter to cloud search
        // let cloud_pvalue = (-profile.forward_lambda as f64
        //     * (cloud_score_bits as f64 - profile.forward_tau as f64))
        //     .exp();

        // if cloud_pvalue >= data.cloud_pvalue_threshold {
        //     continue;
        // }

        self.forward_matrix
            .reuse(target.length, profile.length, bounds);
        self.backward_matrix
            .reuse(target.length, profile.length, bounds);
        self.posterior_matrix
            .reuse(target.length, profile.length, bounds);
        self.optimal_matrix
            .reuse(target.length, profile.length, bounds);

        // we use the forward score to compute the final bit score (later)
        let forward_score_nats = forward(profile, target, &mut self.forward_matrix, bounds);

        let forward_pvalue = (-profile.forward_lambda as f64
            * ((forward_score_nats / std::f32::consts::LN_2) as f64 - profile.forward_tau as f64))
            .exp();

        // TODO: make this work again
        // if forward_pvalue >= data.forward_pvalue_threshold {
        //     continue;
        // }

        backward(profile, target, &mut self.backward_matrix, bounds);

        posterior(
            profile,
            &self.forward_matrix,
            &self.backward_matrix,
            &mut self.posterior_matrix,
            bounds,
        );

        let null1 = length_bias_score(target.length);
        let null2 = composition_bias_score(&self.posterior_matrix, profile, target, bounds);

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

        // TODO: make this work again after I refactor Alignment struct
        // let mut alignment = Alignment::from_trace(&trace, profile, target, &data.score_params);

        // if alignment.evalue <= data.evalue_threshold {
        //     alignment.cell_fraction =
        //         Some(bounds.num_cells() as f32 / (target.length * profile.length) as f32);

        //     data.output.table_format.update_widths(&alignment);
        //     alignments.push(alignment);
        // }
    }
}

struct Pipeline {
    seed: SeedStep,
    cloud_search: CloudSearchStep,
    align: AlignmentStep,
}

impl Pipeline {
    fn run(&mut self, profile: &Profile, target: &Sequence) {
        let seed = self.seed.run(profile, target).unwrap();
        let bounds = self.cloud_search.run(profile, target, &seed).unwrap();
        self.align.run(profile, target, bounds);
    }
}
