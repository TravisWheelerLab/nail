mod cloud_step;
pub use cloud_step::*;

mod seed_step;
pub use seed_step::*;

mod align_step;
pub use align_step::*;

mod output_step;
pub use output_step::*;

use std::collections::HashMap;
use std::sync::{Arc, Mutex};
use std::time::Instant;

use rayon::iter::{IntoParallelRefIterator, IntoParallelRefMutIterator, ParallelIterator};
use thiserror::Error;

use libnail::align::structs::Alignment;
use libnail::structs::{Hmm, Profile, Sequence};

use crate::stats::{FilterStage, Stats, ThreadTimed};

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

#[derive(Default, Clone)]
pub struct PipelineConfig {
    pub cloud_threshold: f64,
    pub forward_threshold: f64,
    pub e_value_threshold: f64,
}

#[derive(Clone)]
pub struct Pipeline {
    pub config: PipelineConfig,
    pub seed: Box<dyn SeedStep + Send + Sync>,
    pub cloud_search: Box<dyn CloudSearchStep + Send + Sync>,
    pub align: Box<dyn AlignStep + Send + Sync>,
    pub output: Arc<Mutex<OutputStep>>,
    pub stats: Stats,
}

impl Pipeline {
    fn run(
        &mut self,
        profile: &mut Profile,
        targets: &HashMap<String, Sequence>,
    ) -> anyhow::Result<Vec<Alignment>> {
        let seeds = self.seed.run(profile, targets);

        let alignments: Vec<Alignment> = match seeds {
            None => vec![],
            Some(seeds) => seeds
                .iter()
                .filter_map(|(target_name, seed)| {
                    self.stats.increment_passed(FilterStage::Seed);

                    let target = targets.get(target_name)?;

                    let now = Instant::now();
                    let maybe_bounds = self.cloud_search.run(profile, target, seed);
                    self.stats
                        .add_threaded_time(ThreadTimed::CloudSearch, now.elapsed());

                    maybe_bounds.map(|bounds| {
                        self.stats.increment_passed(FilterStage::Cloud);
                        self.align.run(profile, target, bounds).unwrap()
                    })
                })
                .collect(),
        };

        alignments.iter().for_each(|ali| {
            self.stats
                .add_threaded_time(ThreadTimed::MemoryInit, ali.times.init);
            self.stats
                .add_threaded_time(ThreadTimed::Forward, ali.times.forward);
        });

        let passed_forward: Vec<_> = alignments
            .into_iter()
            .filter(|ali| ali.boundaries.is_some())
            .collect();

        passed_forward.iter().for_each(|ali| {
            self.stats.increment_passed(FilterStage::Forward);
            self.stats
                .add_threaded_time(ThreadTimed::Backward, ali.times.backward);
            self.stats
                .add_threaded_time(ThreadTimed::Posterior, ali.times.posterior);
            self.stats
                .add_threaded_time(ThreadTimed::Traceback, ali.times.traceback);
            self.stats
                .add_threaded_time(ThreadTimed::NullTwo, ali.times.null_two);
        });

        let mut passed_e_value: Vec<_> = passed_forward
            .into_iter()
            .filter(|ali| ali.scores.e_value <= self.config.e_value_threshold)
            .collect();

        let now = Instant::now();
        match self.output.lock() {
            Ok(mut guard) => guard.write(&mut passed_e_value).unwrap(),
            Err(_) => panic!(),
        };
        self.stats
            .add_threaded_time(ThreadTimed::Output, now.elapsed());

        Ok(passed_e_value)
    }
}

pub fn run_pipeline_profile_to_sequence(
    queries: &mut [Profile],
    targets: &HashMap<String, Sequence>,
    pipeline: &mut Pipeline,
) {
    queries.par_iter_mut().panic_fuse().for_each_with(
        (pipeline.clone(), targets),
        |(pipeline, targets), profile| {
            let now = Instant::now();
            let _ = pipeline.run(profile, targets);
            pipeline
                .stats
                .add_threaded_time(ThreadTimed::Total, now.elapsed())
        },
    );
}

pub fn run_pipeline_sequence_to_sequence(
    queries: &[Sequence],
    targets: &HashMap<String, Sequence>,
    pipeline: &mut Pipeline,
) {
    queries.par_iter().panic_fuse().for_each_with(
        (pipeline.clone(), targets),
        |(pipeline, targets), sequence| {
            let now = Instant::now();
            let mut profile = Hmm::from_blosum_62_and_sequence(sequence)
                .map(|h| Profile::new(&h))
                .expect("failed to build profile from sequence");
            profile.calibrate_tau(200, 100, 0.04);

            pipeline
                .stats
                .add_threaded_time(ThreadTimed::HmmBuild, now.elapsed());

            let _ = pipeline.run(&mut profile, targets);

            pipeline
                .stats
                .add_threaded_time(ThreadTimed::Total, now.elapsed())
        },
    );
}
