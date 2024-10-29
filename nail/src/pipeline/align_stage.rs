use std::time::{Duration, Instant};

use derive_builder::Builder;
use libnail::{
    align::{
        backward, forward, null_one_score, null_two_score, optimal_accuracy, p_value, posterior,
        structs::{Alignment, AlignmentBuilder, DpMatrixSparse, RowBounds, Trace},
        traceback, Bits,
    },
    structs::{Profile, Sequence},
};

use crate::args::SearchArgs;

use super::StageResult;

pub type AlignStageResult = StageResult<Alignment, AlignStageStats>;

impl AlignStageResult {
    pub fn tab_string(&self) -> String {
        match self {
            StageResult::Filtered { stats } => {
                format!(
                    "F {:.2}b {:.1e} {} {}",
                    stats.score.value(),
                    stats.p_value,
                    stats.forward_cells,
                    stats.forward_time.as_nanos(),
                )
            }
            StageResult::Passed { stats, data: ali } => {
                format!(
                    "P {:.2}b {:.1e} {:.1e} {} {} {} {} {} {}",
                    stats.score.value(),
                    stats.p_value,
                    ali.scores.e_value,
                    stats.forward_cells,
                    stats.forward_time.as_nanos(),
                    stats.backward_time.as_nanos(),
                    stats.posterior_time.as_nanos(),
                    stats.optimal_accuracy_time.as_nanos(),
                    stats.null_two_time.as_nanos(),
                )
            }
        }
    }
}

#[derive(Builder, Default)]
#[builder(setter(strip_option), default)]
pub struct AlignStageStats {
    pub forward_cells: usize,
    pub backward_cells: usize,
    pub score: Bits,
    pub p_value: f64,
    pub memory_init_time: Duration,
    pub forward_time: Duration,
    pub backward_time: Duration,
    pub posterior_time: Duration,
    pub optimal_accuracy_time: Duration,
    pub traceback_time: Duration,
    pub null_two_time: Duration,
}

impl AlignStageStatsBuilder {
    fn add_memory_init_time(&mut self, duration: Duration) {
        match self.memory_init_time {
            Some(ref mut time) => *time += duration,
            None => {
                self.memory_init_time(duration);
            }
        }
    }
}

#[derive(Clone)]
pub struct AlignConfig {
    pub do_null_two: bool,
}

impl Default for AlignConfig {
    fn default() -> Self {
        Self { do_null_two: true }
    }
}

pub trait AlignStage: dyn_clone::DynClone {
    fn run(
        &mut self,
        profile: &mut Profile,
        target: &Sequence,
        bounds: &RowBounds,
    ) -> StageResult<Alignment, AlignStageStats>;
}

dyn_clone::clone_trait_object!(AlignStage);

#[derive(Default, Clone)]
pub struct DefaultAlignStage {
    forward_matrix: DpMatrixSparse,
    backward_matrix: DpMatrixSparse,
    posterior_matrix: DpMatrixSparse,
    optimal_matrix: DpMatrixSparse,
    forward_p_value_threshold: f64,
    target_count: usize,
    config: AlignConfig,
}

impl DefaultAlignStage {
    pub fn new(args: &SearchArgs) -> Self {
        Self {
            target_count: match args.nail_args.target_database_size {
                Some(size) => size,
                None => panic!(),
            },
            forward_p_value_threshold: args.nail_args.forward_pvalue_threshold,
            ..Default::default()
        }
    }
}

impl AlignStage for DefaultAlignStage {
    fn run(
        &mut self,
        profile: &mut Profile,
        target: &Sequence,
        bounds: &RowBounds,
    ) -> StageResult<Alignment, AlignStageStats> {
        let mut stats = AlignStageStatsBuilder::default();

        // configuring for the target length adjusts special state transitions
        profile.configure_for_target_length(target.length);

        let now = Instant::now();
        self.forward_matrix
            .reuse(target.length, profile.length, bounds);
        stats.memory_init_time(now.elapsed());

        // we use the forward score to compute the final bit score (later)
        let now = Instant::now();
        let forward_score = forward(profile, target, &mut self.forward_matrix, bounds).to_bits()
            // the denominator is the null one score
            - null_one_score(target.length);
        stats.forward_time(now.elapsed());
        stats.forward_cells(bounds.num_cells);

        // for now we compute the P-value for filtering purposes
        let forward_p_value = p_value(forward_score, profile.forward_lambda, profile.forward_tau);
        stats.score(forward_score);
        stats.p_value(forward_p_value);

        if forward_p_value >= self.forward_p_value_threshold {
            return StageResult::Filtered {
                stats: stats.build().unwrap(),
            };
        }

        let now = Instant::now();
        self.backward_matrix
            .reuse(target.length, profile.length, bounds);
        self.posterior_matrix
            .reuse(target.length, profile.length, bounds);
        self.optimal_matrix
            .reuse(target.length, profile.length, bounds);
        stats.add_memory_init_time(now.elapsed());

        let now = Instant::now();
        backward(profile, target, &mut self.backward_matrix, bounds);
        stats.backward_time(now.elapsed());
        stats.backward_cells(bounds.num_cells);

        let now = Instant::now();
        posterior(
            profile,
            &self.forward_matrix,
            &self.backward_matrix,
            &mut self.posterior_matrix,
            bounds,
        );
        stats.posterior_time(now.elapsed());

        let now = Instant::now();
        optimal_accuracy(
            profile,
            &self.posterior_matrix,
            &mut self.optimal_matrix,
            bounds,
        );
        stats.optimal_accuracy_time(now.elapsed());

        let now = Instant::now();
        let mut trace = Trace::new(target.length, profile.length);
        traceback(
            profile,
            &self.posterior_matrix,
            &self.optimal_matrix,
            &mut trace,
            bounds.target_end,
        );
        stats.traceback_time(now.elapsed());

        let null_two_score = if self.config.do_null_two {
            let now = Instant::now();
            let score = Some(null_two_score(
                &self.posterior_matrix,
                profile,
                target,
                bounds,
            ));
            stats.null_two_time(now.elapsed());
            score
        } else {
            None
        };

        StageResult::Passed {
            data: AlignmentBuilder::default()
                .with_profile(profile)
                .with_target(target)
                .with_database_size(self.target_count)
                .with_cell_count(bounds.num_cells)
                .with_forward_score(forward_score)
                .with_trace(&trace)
                .with_null_two(null_two_score)
                .build()
                .unwrap(),
            stats: stats.build().unwrap(),
        }
    }
}
