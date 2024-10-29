use std::time::{Duration, Instant};

use derive_builder::Builder;
use libnail::{
    align::{
        cloud_score, cloud_search_backward, cloud_search_forward, p_value,
        structs::{AntiDiagonalBounds, CloudMatrixLinear, RowBounds, Seed},
        CloudSearchParams, Nats,
    },
    structs::{Profile, Sequence},
};

use crate::args::SearchArgs;

use super::StageResult;

pub type CloudStageResult = StageResult<RowBounds, CloudStageStats>;

impl CloudStageResult {
    pub fn tab_string(&self) -> String {
        match self {
            StageResult::Filtered { stats } => {
                format!(
                    "F {:.2}n {:.1e} {} {} {} {}",
                    stats.score.value(),
                    stats.p_value,
                    stats.forward_cells,
                    stats.backward_cells,
                    stats.forward_time.as_nanos(),
                    stats.backward_time.as_nanos(),
                )
            }
            StageResult::Passed {
                stats,
                data: bounds,
            } => {
                format!(
                    "P {:.2}b {:.1e} {} {} {} {} {}",
                    stats.score.value(),
                    stats.p_value,
                    stats.forward_cells,
                    stats.backward_cells,
                    bounds.num_cells,
                    stats.forward_time.as_nanos(),
                    stats.backward_time.as_nanos(),
                )
            }
        }
    }
}

#[derive(Builder, Default)]
#[builder(setter(strip_option), default)]
pub struct CloudStageStats {
    pub forward_cells: usize,
    pub backward_cells: usize,
    pub score: Nats,
    pub p_value: f64,
    pub memory_init_time: Duration,
    pub forward_time: Duration,
    pub backward_time: Duration,
    pub merge_time: Duration,
    pub trim_time: Duration,
    pub reorient_time: Duration,
}

pub trait CloudSearchStage: dyn_clone::DynClone {
    fn run(&mut self, profile: &Profile, target: &Sequence, seed: &Seed) -> CloudStageResult;
}

dyn_clone::clone_trait_object!(CloudSearchStage);

#[derive(Default, Clone)]
pub struct FullDpCloudSearchStage {}

impl CloudSearchStage for FullDpCloudSearchStage {
    fn run(&mut self, profile: &Profile, target: &Sequence, _seed: &Seed) -> CloudStageResult {
        let mut row_bounds = RowBounds::default();
        row_bounds.fill_rectangle(1, 1, target.length, profile.length);

        StageResult::Passed {
            data: row_bounds,
            stats: CloudStageStats::default(),
        }
    }
}

#[derive(Default, Clone)]
pub struct DefaultCloudSearchStage {
    cloud_matrix: CloudMatrixLinear,
    forward_bounds: AntiDiagonalBounds,
    reverse_bounds: AntiDiagonalBounds,
    params: CloudSearchParams,
    p_value_threshold: f64,
}

impl DefaultCloudSearchStage {
    pub fn new(args: &SearchArgs) -> Self {
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

impl CloudSearchStage for DefaultCloudSearchStage {
    fn run(&mut self, profile: &Profile, target: &Sequence, seed: &Seed) -> CloudStageResult {
        let mut stats = CloudStageStatsBuilder::default();

        let now = Instant::now();
        self.cloud_matrix.reuse(profile.length);
        self.forward_bounds.reuse(target.length, profile.length);
        self.reverse_bounds.reuse(target.length, profile.length);
        let mut row_bounds = RowBounds::new(target.length);
        stats.memory_init_time(now.elapsed());

        let now = Instant::now();
        let forward_results = cloud_search_forward(
            profile,
            target,
            seed,
            &mut self.cloud_matrix,
            &self.params,
            &mut self.forward_bounds,
        );
        stats.forward_time(now.elapsed());
        stats.forward_cells(forward_results.num_cells_computed);

        let now = Instant::now();
        let backward_results = cloud_search_backward(
            profile,
            target,
            seed,
            &mut self.cloud_matrix,
            &self.params,
            &mut self.reverse_bounds,
        );
        stats.backward_time(now.elapsed());
        stats.backward_cells(backward_results.num_cells_computed);

        let cloud_score = cloud_score(&forward_results, &backward_results);
        let cloud_p_value = p_value(cloud_score, profile.forward_lambda, profile.forward_tau);

        stats.score(cloud_score);
        stats.p_value(cloud_p_value);

        if cloud_p_value >= self.p_value_threshold {
            return StageResult::Filtered {
                stats: stats.build().unwrap(),
            };
        }

        let now = Instant::now();
        self.forward_bounds.merge(&self.reverse_bounds);
        stats.merge_time(now.elapsed());

        self.forward_bounds.square_corners();

        let now = Instant::now();
        let trim_result = self.forward_bounds.trim_wings();
        stats.trim_time(now.elapsed());

        match trim_result {
            Ok(_) => {
                let now = Instant::now();
                row_bounds.fill_from_anti_diagonal_bounds(&self.forward_bounds);
                stats.reorient_time(now.elapsed());
            }
            // TODO: probably want to do something else/extra here
            Err(_) => {
                row_bounds.fill_rectangle(
                    seed.target_start,
                    seed.profile_start,
                    seed.target_end,
                    seed.profile_end,
                );
            }
        }

        StageResult::Passed {
            data: row_bounds,
            stats: stats.build().unwrap(),
        }
    }
}
