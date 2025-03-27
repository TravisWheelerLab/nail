use std::time::{Duration, Instant};

use derive_builder::Builder;
use libnail::{
    align::{
        cloud_score, cloud_search_backward, cloud_search_forward, cloud_search_test, forward_test,
        null_one_score, p_value,
        structs::{
            AdMatrixSparse, Cloud, CloudMatrixLinear, CloudMatrixQuadratic, DpMatrixFlat,
            RowBounds, Seed,
        },
        CloudSearchParams, Nats,
    },
    structs::{sequence, Profile, Sequence},
};

use crate::args::SearchArgs;

use super::StageResult;

pub type CloudStageResult = StageResult<RowBounds, CloudStageStats>;

impl CloudStageResult {
    pub fn tab_string(&self) -> String {
        let (stats, pass_str) = match self {
            StageResult::Filtered { stats } => (stats, "F"),
            StageResult::Passed { stats, .. } => (stats, "P"),
        };
        format!(
            "{} {:.2}b {:.1e} {} {} {} {}",
            pass_str,
            stats.score.to_bits().value(),
            stats.p_value,
            stats.forward_cells,
            stats.backward_cells,
            stats.forward_time.as_nanos(),
            stats.backward_time.as_nanos(),
        )
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

pub trait CloudSearchStage: dyn_clone::DynClone + Send + Sync {
    fn run(&mut self, profile: &Profile, target: &Sequence, seed: &Seed) -> CloudStageResult;
}

dyn_clone::clone_trait_object!(CloudSearchStage);

#[derive(Default, Clone)]
pub struct TmpDebugCloudSearchStage {}

impl CloudSearchStage for TmpDebugCloudSearchStage {
    fn run(&mut self, profile: &Profile, target: &Sequence, _: &Seed) -> CloudStageResult {
        let mut profile = profile.clone();
        profile.configure_for_target_length(target.length);

        // let mut fwd_lin = CloudMatrixLinear::new(profile.length);
        // let mut bwd_lin = CloudMatrixLinear::new(profile.length);
        // cloud_search_test(profile, target, seed, &mut fwd_lin, &mut bwd_lin);

        // let mut fwd_cloud_quad = CloudMatrixQuadratic::new(target.length, profile.length);
        // let mut bwd_cloud_quad = CloudMatrixQuadratic::new(target.length, profile.length);

        // cloud_search_test(
        //     &mut profile,
        //     target,
        //     &mut fwd_cloud_quad,
        //     &mut bwd_cloud_quad,
        // );

        let mut fwd_matrix_og = DpMatrixFlat::new(target.length, profile.length);
        let mut bwd_matrix_og = DpMatrixFlat::new(target.length, profile.length);
        forward_test(
            &mut profile,
            target,
            &mut fwd_matrix_og,
            &mut bwd_matrix_og,
            "1",
        );

        let mut bounds = Cloud::new(target.length, profile.length);
        bounds.fill_rectangle(0, 0, target.length, profile.length);
        let mut fwd_matrix_new = AdMatrixSparse::from_cloud(&bounds);
        let mut bwd_matrix_new = AdMatrixSparse::from_cloud(&bounds);

        forward_test(
            &mut profile,
            target,
            &mut fwd_matrix_new,
            &mut bwd_matrix_new,
            "2",
        );

        fwd_matrix_new.layout();

        StageResult::Filtered {
            stats: CloudStageStatsBuilder::default().build().unwrap(),
        }
    }
}

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
    forward_bounds: Cloud,
    reverse_bounds: Cloud,
    params: CloudSearchParams,
    p_value_threshold: f64,
}

impl DefaultCloudSearchStage {
    pub fn new(args: &SearchArgs) -> Self {
        Self {
            params: CloudSearchParams {
                gamma: args.pipeline_args.gamma,
                alpha: args.pipeline_args.alpha,
                beta: args.pipeline_args.beta,
            },
            p_value_threshold: args.pipeline_args.cloud_pvalue_threshold,
            ..Default::default()
        }
    }
}

impl CloudSearchStage for DefaultCloudSearchStage {
    fn run(&mut self, profile: &Profile, target: &Sequence, seed: &Seed) -> CloudStageResult {
        let now = Instant::now();
        let mut stats = CloudStageStatsBuilder::default();
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

        // cloud search does not compute any special state scores
        let raw_cloud_score = cloud_score(&forward_results, &backward_results);

        // TODO: put this score adjustment somewhere else
        let target_length = target.length as f32;
        let profile_length = profile.length as f32;
        let cloud_score_adjustment = Nats(
            // L * log(L / (L + 2))          [ N->N | C->C ] non-homologous states
            target_length * (target_length / (target_length + 2.0)).ln()
                // NOTE: no longer adding this here
                // // 2 * log(2 / (L + 2))      [ N->B + C->T ] start/end transitions
                // + 2.0 * (2.0 / (target_length + 2.0)).ln()
                // log(2 / (M * (M + 1)))    [ B->M + M->E ] core model entry approximation
                + (2.0 / (profile_length * (profile_length + 1.0))).ln(),
        );

        // the null one is the denominator in the probability ratio
        let null_one = null_one_score(target.length);

        let cloud_score = raw_cloud_score + cloud_score_adjustment - null_one;

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
