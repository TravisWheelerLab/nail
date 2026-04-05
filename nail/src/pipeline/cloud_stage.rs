use std::time::{Duration, Instant};

use derive_builder::Builder;
use libnail::{
    align::{
        cloud_search_bwd, cloud_search_fwd, null_one_score, p_value,
        structs::{
            AdMatrixLinear,
            BackgroundState::{C as BgC, N as BgN},
            Cloud,
            Relationship::{Disjoint, Intersecting},
            RowBounds, Seed,
        },
        Bits, CloudSearchParams, Nats,
    },
    structs::{Profile, Sequence},
};

use crate::args::SearchArgs;

use super::StageResult;

pub type CloudStageResult = StageResult<RowBounds, CloudStageStats>;

#[derive(Builder, Default)]
#[builder(setter(strip_option), default)]
pub struct CloudStageStats {
    pub forward_cells: usize,
    pub backward_cells: usize,
    pub score: Bits,
    pub p_value: f64,
    pub memory_init_time: Duration,
    pub forward_time: Duration,
    pub backward_time: Duration,
    pub merge_time: Duration,
    pub trim_time: Duration,
    pub reorient_time: Duration,
}

pub trait CloudSearchStage: dyn_clone::DynClone + Send + Sync {
    fn run(&mut self, prf: &Profile, seq: &Sequence, seed: &Seed) -> CloudStageResult;
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
    mx: AdMatrixLinear,
    fwd_cloud: Cloud,
    bwd_cloud: Cloud,
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
    fn run(&mut self, prf: &Profile, seq: &Sequence, seed: &Seed) -> CloudStageResult {
        let now = Instant::now();
        let mut stats = CloudStageStatsBuilder::default();
        let mut row_bounds = RowBounds::new(seq.length);
        stats.memory_init_time(now.elapsed());

        let mut fwd_results = None;
        let mut bwd_results = None;
        let mut raw_cloud_score = None;
        let mut raw_bwd_cloud_score = None;
        let mut clouds_intersected = false;

        let null_one = null_one_score(seq.length);

        for attempts in 0..5 {
            let now = Instant::now();
            self.mx.reuse(seq.length);
            self.fwd_cloud.reuse(seq.length, prf.length);
            self.bwd_cloud.reuse(seq.length, prf.length);
            stats.memory_init_time(now.elapsed());

            let params = self.params.scale(1.0 + (attempts as f32 * 0.5));
            let now = Instant::now();
            fwd_results = Some(cloud_search_fwd(
                prf,
                seq,
                seed,
                &mut self.mx,
                &params,
                &mut self.fwd_cloud,
            ));
            stats.forward_time(now.elapsed());

            raw_cloud_score = Some(Nats(
                self.mx[(BgC, seed.seq_end)]
                    + prf.special_transition_score(Profile::C_IDX, Profile::SPECIAL_MOVE_IDX),
            ));

            // DEBUG: print seed coordinates and C-state scores
            if attempts == 0 {
                eprintln!(
                    "DEBUG seed: seq={}..{} prf={}..{} seq_len={}",
                    seed.seq_start, seed.seq_end, seed.prf_start, seed.prf_end, seq.length
                );
                for pos in seed.seq_start..=seed.seq_end {
                    eprintln!(
                        "DEBUG C[{}] = {:.4}",
                        pos,
                        self.mx[(BgC, pos)]
                    );
                }
            }

            let now = Instant::now();
            self.mx.reuse(seq.length);
            stats.memory_init_time(now.elapsed());

            let now = Instant::now();
            bwd_results = Some(cloud_search_bwd(
                prf,
                seq,
                seed,
                &mut self.mx,
                &params,
                &mut self.bwd_cloud,
            ));
            stats.backward_time(now.elapsed());

            raw_bwd_cloud_score = Some(Nats(self.mx[(BgN, seed.seq_start)]));

            // DEBUG: print backward N-state score
            if attempts == 0 {
                eprintln!(
                    "DEBUG N[{}] = {:.4}",
                    seed.seq_start,
                    self.mx[(BgN, seed.seq_start)]
                );
            }

            match self.fwd_cloud.anti_diagonal_relationship(&self.bwd_cloud) {
                Disjoint(_) => {
                    if attempts >= 4 {
                        // clouds never intersected; check backward score before giving up
                        let bwd_score =
                            (raw_bwd_cloud_score.unwrap() - null_one).to_bits();
                        let bwd_p_value =
                            p_value(bwd_score, prf.fwd_lambda, prf.fwd_tau);
                        if bwd_p_value < self.p_value_threshold {
                            // backward pass found real signal; proceed with seed-based bounds
                            break;
                        }
                        return StageResult::Filtered {
                            stats: stats.build().unwrap(),
                        };
                    }
                }
                Intersecting(_) => {
                    clouds_intersected = true;
                    break;
                }
            };
        }

        let (fwd_results, bwd_results, raw_cloud_score, raw_bwd_cloud_score) =
            match (fwd_results, bwd_results, raw_cloud_score, raw_bwd_cloud_score) {
                (Some(f), Some(b), Some(s), Some(t)) => (f, b, s, t),
                _ => unreachable!("finished cloud search without results"),
            };

        stats.backward_cells(bwd_results.num_cells_computed);
        stats.forward_cells(fwd_results.num_cells_computed);

        // use the better of forward C-state or backward N-state score
        let best_raw_score = if raw_bwd_cloud_score > raw_cloud_score {
            raw_bwd_cloud_score
        } else {
            raw_cloud_score
        };

        let cloud_score = (best_raw_score - null_one).to_bits();

        let cloud_p_value = p_value(cloud_score, prf.fwd_lambda, prf.fwd_tau);
        stats.score(cloud_score);
        stats.p_value(cloud_p_value);

        if cloud_p_value >= self.p_value_threshold {
            return StageResult::Filtered {
                stats: stats.build().unwrap(),
            };
        }

        if clouds_intersected {
            let now = Instant::now();
            self.fwd_cloud.merge(&self.bwd_cloud);
            stats.merge_time(now.elapsed());

            self.fwd_cloud.square_corners();

            let now = Instant::now();
            let trim_result = self.fwd_cloud.trim_wings();
            stats.trim_time(now.elapsed());

            match trim_result {
                Ok(_) => {
                    let now = Instant::now();
                    row_bounds.fill_from_cloud(&self.fwd_cloud);
                    stats.reorient_time(now.elapsed());
                }
                // TODO: probably want to do something else/extra here
                Err(_) => {
                    row_bounds.fill_rectangle(
                        seed.seq_start,
                        seed.prf_start,
                        seed.seq_end,
                        seed.prf_end,
                    );
                }
            }
        } else {
            // clouds were disjoint but backward score passed — use seed bounds directly
            row_bounds.fill_rectangle(
                seed.seq_start,
                seed.prf_start,
                seed.seq_end,
                seed.prf_end,
            );
        }

        StageResult::Passed {
            data: row_bounds,
            stats: stats.build().unwrap(),
        }
    }
}
