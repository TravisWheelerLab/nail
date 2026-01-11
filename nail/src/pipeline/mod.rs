mod cloud_stage;
pub use cloud_stage::*;

mod seed_stage;
pub use seed_stage::*;

mod align_stage;
pub use align_stage::*;

mod output_stage;
pub use output_stage::*;

use std::cell::RefCell;
use std::time::Instant;

use rayon::iter::{IntoParallelRefIterator, ParallelIterator};
use thiserror::Error;
use thread_local::ThreadLocal;

use libnail::{
    align::{structs::Seed, Bits},
    structs::Profile,
};

use crate::{
    io::{Fasta, P7Hmm},
    stats::{Stats, ThreadedTimed},
};

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
    pub seed_result: SeedStageResult,
    pub cloud_result: Option<CloudStageResult>,
    pub align_result: Option<AlignStageResult>,
}

pub trait TableDisplay {
    fn cell(&self) -> String;
}

impl<T> TableDisplay for Option<T>
where
    T: TableDisplay,
{
    fn cell(&self) -> String {
        match self {
            Some(v) => v.cell(),
            None => "-".to_string(),
        }
    }
}

impl TableDisplay for f64 {
    fn cell(&self) -> String {
        format!("{:.1e}", self)
    }
}

impl TableDisplay for Bits {
    fn cell(&self) -> String {
        format!("{:.1}", self.0)
    }
}

impl PipelineResult {
    pub fn stat_string(&self) -> String {
        let seed_stats = self.seed_result.stats();
        let (seed_s, seed_e) = (Some(seed_stats.score), Some(seed_stats.e_value));

        let (cloud_s, cloud_p) = self
            .cloud_result
            .as_ref()
            .map(|r| (Some(r.stats().score), Some(r.stats().p_value)))
            .unwrap_or_default();

        let (align_s, align_p) = self
            .align_result
            .as_ref()
            .map(|r| (Some(r.stats().score), Some(r.stats().p_value)))
            .unwrap_or_default();

        format!(
            "{} {} {} {} {} {} {} {} {} {}",
            self.profile_name,
            self.target_name,
            self.profile_length,
            self.target_length,
            seed_s.cell(),
            seed_e.cell(),
            cloud_s.cell(),
            cloud_p.cell(),
            align_s.cell(),
            align_p.cell(),
        )
    }
}

#[derive(Clone)]
pub struct Pipeline {
    pub targets: Fasta,
    pub seed: Box<dyn SeedStage>,
    pub cloud_search: Box<dyn CloudSearchStage>,
    pub align: Box<dyn AlignStage>,
    pub output: OutputStage,
    pub stats: Stats,
}

impl Pipeline {
    fn run(&mut self, profile: &mut Profile) -> anyhow::Result<()> {
        let seeds = self.seed.run(profile);

        let pipeline_results: Vec<PipelineResult> = match seeds {
            None => return Ok(()),
            Some(seeds) => seeds
                .into_iter()
                .filter_map(|(seq_name, seed)| {
                    let target = match self.targets.get(&seq_name) {
                        Some(target) => target,
                        // TODO: probably return an error here instead
                        None => return None,
                    };

                    // configuring for the target length adjusts special state transitions
                    profile.configure_for_target_length(target.length);

                    let cloud_result = self.cloud_search.run(profile, &target, &seed);

                    let align_result = match cloud_result {
                        StageResult::Passed {
                            data: ref bounds, ..
                        } => Some(self.align.run(profile, &target, bounds)),
                        StageResult::Filtered { .. } => None,
                    };

                    Some(PipelineResult {
                        profile_name: profile.name.clone(),
                        target_name: target.name.clone(),
                        profile_length: profile.length,
                        target_length: target.length,
                        seed_result: StageResult::Passed {
                            data: seed.clone(),
                            stats: SeedStageStats {
                                score: Bits(seed.score),
                                e_value: seed.e_value,
                            },
                        },
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

    fn run2(
        &mut self,
        profile: &mut Profile,
        seq_name: &str,
        seed: &Seed,
    ) -> anyhow::Result<PipelineResult> {
        let target = match self.targets.get(seq_name) {
            Some(target) => target,
            // TODO: probably return an error here instead
            None => panic!(),
        };

        // configuring for the target length adjusts special state transitions
        profile.configure_for_target_length(target.length);

        let cloud_result = self.cloud_search.run(profile, &target, seed);

        let align_result = match cloud_result {
            StageResult::Passed {
                data: ref bounds, ..
            } => Some(self.align.run(profile, &target, bounds)),
            StageResult::Filtered { .. } => None,
        };

        Ok(PipelineResult {
            profile_name: profile.name.clone(),
            target_name: target.name.clone(),
            profile_length: profile.length,
            target_length: target.length,
            seed_result: StageResult::Passed {
                data: seed.clone(),
                stats: SeedStageStats {
                    score: Bits(seed.score),
                    e_value: seed.e_value,
                },
            },
            cloud_result: Some(cloud_result),
            align_result,
        })
    }
}

pub fn run_pipeline_profile_to_sequence(profiles: &mut P7Hmm, pipeline: &mut Pipeline) {
    let tl_pipeline: ThreadLocal<RefCell<Pipeline>> = ThreadLocal::new();

    let prf_names = profiles
        .names_iter()
        .map(|s| s.to_string())
        .collect::<Vec<_>>();

    for n in prf_names.iter() {
        let prf = profiles.get(n).expect("failed to retreive profile");
        if let Some(seeds) = pipeline.seed.run(&prf) {
            let tl_prf: ThreadLocal<RefCell<Profile>> = ThreadLocal::new();

            let res = seeds
                .par_iter()
                .panic_fuse()
                .map(|(seq_name, seed)| {
                    let now = Instant::now();

                    let mut prf = tl_prf.get_or(|| RefCell::new(prf.clone())).borrow_mut();
                    let mut pipeline = tl_pipeline
                        .get_or(|| RefCell::new(pipeline.clone()))
                        .borrow_mut();

                    let res = pipeline.run2(&mut prf, seq_name, seed);

                    pipeline
                        .stats
                        .add_threaded_time(ThreadedTimed::Total, now.elapsed());

                    res
                })
                .collect::<Result<Vec<_>, _>>()
                .unwrap();

            let output_stats = pipeline.output.run(&res).unwrap();
            pipeline.stats.add_sample(&res, &output_stats);
        }
    }
}

pub fn run_pipeline_sequence_to_sequence(queries: &Fasta, pipeline: &mut Pipeline) {
    let thread_local_pipeline: ThreadLocal<RefCell<Pipeline>> = ThreadLocal::new();

    queries.par_iter().panic_fuse().for_each(|seq| {
        let now = Instant::now();

        let mut pipeline = thread_local_pipeline
            .get_or(|| RefCell::new(pipeline.clone()))
            .borrow_mut();

        let mut profile =
            Profile::from_blosum_62_and_seq(&seq).expect("failed to build profile from sequence");

        pipeline
            .stats
            .add_threaded_time(ThreadedTimed::HmmBuild, now.elapsed());

        let _ = pipeline.run(&mut profile);

        pipeline
            .stats
            .add_threaded_time(ThreadedTimed::Total, now.elapsed())
    });
}
