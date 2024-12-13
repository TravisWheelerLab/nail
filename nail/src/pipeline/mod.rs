mod cloud_stage;
pub use cloud_stage::*;

mod seed_stage;
use indicatif::ProgressBar;
pub use seed_stage::*;

mod align_stage;
pub use align_stage::*;

mod output_stage;
pub use output_stage::*;
use thread_local::ThreadLocal;

use std::cell::RefCell;
use std::collections::HashMap;
use std::sync::Arc;
use std::time::Instant;

use rayon::iter::{IntoParallelRefIterator, IntoParallelRefMutIterator, ParallelIterator};
use thiserror::Error;

use libnail::structs::{Hmm, Profile, Sequence};

use crate::stats::{Stats, ThreadedTimed};

#[derive(Error, Debug)]
#[error("no profile with name: {profile_name}")]
pub struct ProfileNotFoundError {
    profile_name: String,
}

#[derive(Error, Debug)]
#[error("no target with name: {target_name}")]
pub struct TargetNotFoundError {
    target_name: String,
}

pub enum StageResult<D, S> {
    Filtered { stats: S },
    Passed { data: D, stats: S },
}

impl<D, S> StageResult<D, S> {
    pub fn stats(&self) -> &S {
        match self {
            StageResult::Passed { data: _, stats } => stats,
            StageResult::Filtered { stats } => stats,
        }
    }
}

pub struct PipelineResult {
    pub profile_name: String,
    pub target_name: String,
    pub profile_length: usize,
    pub target_length: usize,
    pub cloud_result: Option<CloudStageResult>,
    pub align_result: Option<AlignStageResult>,
}

impl PipelineResult {
    pub fn tab_string(&self) -> String {
        format!(
            "({} {} {} {}) ({}) ({})",
            self.profile_name,
            self.target_name,
            self.profile_length,
            self.target_length,
            match self.cloud_result.as_ref() {
                Some(result) => result.tab_string(),
                None => "".to_string(),
            },
            match self.align_result.as_ref() {
                Some(result) => result.tab_string(),
                None => "".to_string(),
            },
        )
    }
}

#[derive(Clone)]
pub struct Pipeline {
    pub seed: Box<dyn SeedStage>,
    pub cloud_search: Box<dyn CloudSearchStage>,
    pub align: Box<dyn AlignStage>,
    pub output: OutputStage,
    pub stats: Stats,
}

impl Pipeline {
    fn run(
        &mut self,
        profile: &mut Profile,
        targets: &HashMap<String, Sequence>,
    ) -> anyhow::Result<()> {
        let seeds = self.seed.run(profile, targets);

        let pipeline_results: Vec<PipelineResult> = match seeds {
            None => return Ok(()),
            Some(seeds) => seeds
                .iter()
                .filter_map(|(target_name, seed)| {
                    let target = match targets.get(target_name) {
                        Some(target) => target,
                        // TODO: probably return an error here instead
                        None => return None,
                    };

                    let cloud_result = self.cloud_search.run(profile, target, seed);

                    let align_result = match cloud_result {
                        StageResult::Passed {
                            data: ref bounds, ..
                        } => Some(self.align.run(profile, target, bounds)),
                        StageResult::Filtered { .. } => None,
                    };

                    Some(PipelineResult {
                        profile_name: profile.name.clone(),
                        target_name: target.name.clone(),
                        profile_length: profile.length,
                        target_length: target.length,
                        cloud_result: Some(cloud_result),
                        align_result,
                    })
                })
                .collect(),
        };

        let output_stats = self.output.run(&pipeline_results)?;
        self.stats.add_sample(&pipeline_results, &output_stats);

        Ok(())
    }
}

pub fn run_pipeline_profile_to_sequence(
    queries: &mut [Profile],
    targets: &HashMap<String, Sequence>,
    pipeline: &mut Pipeline,
) {
    let bar = Arc::new(ProgressBar::new(queries.len() as u64));
    let thread_local_pipeline: ThreadLocal<RefCell<Pipeline>> = ThreadLocal::new();

    queries
        .par_iter_mut()
        .panic_fuse()
        .for_each_with(targets, |targets, profile| {
            let now = Instant::now();
            let mut pipeline = thread_local_pipeline
                .get_or(|| RefCell::new(pipeline.clone()))
                .borrow_mut();

            let _ = pipeline.run(profile, targets);

            bar.inc(1);

            pipeline
                .stats
                .add_threaded_time(ThreadedTimed::Total, now.elapsed())
        });

    bar.finish();
}

pub fn run_pipeline_sequence_to_sequence(
    queries: &[Sequence],
    targets: &HashMap<String, Sequence>,
    pipeline: &mut Pipeline,
) {
    let bar = Arc::new(ProgressBar::new(queries.len() as u64));
    let thread_local_pipeline: ThreadLocal<RefCell<Pipeline>> = ThreadLocal::new();

    queries
        .par_iter()
        .panic_fuse()
        .for_each_with(targets, |targets, sequence| {
            let now = Instant::now();

            let mut pipeline = thread_local_pipeline
                .get_or(|| RefCell::new(pipeline.clone()))
                .borrow_mut();

            let mut profile = Hmm::from_blosum_62_and_sequence(sequence)
                .map(|h| Profile::new(&h))
                .expect("failed to build profile from sequence");
            profile.calibrate_tau(200, 100, 0.04);

            pipeline
                .stats
                .add_threaded_time(ThreadedTimed::HmmBuild, now.elapsed());

            let _ = pipeline.run(&mut profile, targets);

            bar.inc(1);

            pipeline
                .stats
                .add_threaded_time(ThreadedTimed::Total, now.elapsed())
        });

    bar.finish();
}
