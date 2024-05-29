use libnail::{
    align::{
        cloud_score, cloud_search_backward, cloud_search_forward, p_value,
        structs::{AntiDiagonalBounds, CloudMatrixLinear, RowBounds, Seed},
        CloudSearchParams,
    },
    structs::{Profile, Sequence},
};

use crate::args::AlignArgs;

pub trait CloudSearchStep: dyn_clone::DynClone {
    fn run(&mut self, profile: &Profile, target: &Sequence, seed: &Seed) -> Option<&RowBounds>;
}

dyn_clone::clone_trait_object!(CloudSearchStep);

#[derive(Default, Clone)]
pub struct FullDpCloudSearchStep {
    row_bounds: RowBounds,
}

impl CloudSearchStep for FullDpCloudSearchStep {
    fn run(&mut self, profile: &Profile, target: &Sequence, _seed: &Seed) -> Option<&RowBounds> {
        self.row_bounds
            .fill_rectangle(1, 1, target.length, profile.length);
        Some(&self.row_bounds)
    }
}

#[derive(Default, Clone)]
pub struct DefaultCloudSearchStep {
    cloud_matrix: CloudMatrixLinear,
    forward_bounds: AntiDiagonalBounds,
    reverse_bounds: AntiDiagonalBounds,
    row_bounds: RowBounds,
    params: CloudSearchParams,
    p_value_threshold: f64,
}

impl DefaultCloudSearchStep {
    pub fn new(args: &AlignArgs) -> Self {
        Self {
            params: CloudSearchParams {
                gamma: args.nail_args.gamma,
                alpha: args.nail_args.alpha,
                beta: args.nail_args.beta,
            },
            p_value_threshold: args.nail_args.cloud_pvalue_threshold,
            ..Default::default()
        }
    }
}

impl CloudSearchStep for DefaultCloudSearchStep {
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

        let cloud_score = cloud_score(&forward_scores, &reverse_scores);

        let cloud_p_value = p_value(cloud_score, profile.forward_lambda, profile.forward_tau);

        if cloud_p_value >= self.p_value_threshold {
            return None;
        }

        let bounds_intersect =
            self.forward_bounds.max_anti_diagonal_idx >= self.reverse_bounds.min_anti_diagonal_idx;

        // TODO: clean up this mess
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
