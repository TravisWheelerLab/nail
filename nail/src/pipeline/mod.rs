mod cloud_stage;
use anyhow::Context;
pub use cloud_stage::*;

mod seed_stage;
pub use seed_stage::*;

mod align_stage;
pub use align_stage::*;

mod output_stage;
pub use output_stage::*;

use std::collections::HashMap;
use std::time::Instant;
use std::{cell::RefCell, sync::Arc};

use rayon::{iter::ParallelIterator, slice::ParallelSlice};

use thread_local::ThreadLocal;

use libnail::{
    align::{structs::Seed, Bits},
    structs::Profile,
};

use crate::{
    io::{Fasta, Seeds2},
    stats::{Stats, ThreadedTimed},
};

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
    pub profiles: Arc<HashMap<String, Profile>>,
    pub prf: Option<Profile>,
    pub targets: Fasta,
    pub cloud_search: Box<dyn CloudSearchStage>,
    pub align: Box<dyn AlignStage>,
    pub output: OutputStage,
    pub stats: Stats,
}

impl Pipeline {
    fn run(&mut self, seed: &Seed) -> anyhow::Result<PipelineResult> {
        match self.prf {
            Some(ref prf) => {
                if prf.name != seed.prf {
                    self.prf = self.profiles.get(&seed.prf).cloned();
                }
            }
            None => {
                self.prf = self.profiles.get(&seed.prf).cloned();
            }
        }

        let prf = self.prf.as_mut().expect("no prf");
        let seq = self.targets.get(&seed.seq).context("no seq")?;

        // configuring for the target length adjusts special state transitions
        prf.configure_for_target_length(seq.length);

        let cloud_result = self.cloud_search.run(prf, &seq, seed);

        let align_result = match cloud_result {
            StageResult::Passed {
                data: ref bounds, ..
            } => Some(self.align.run(prf, &seq, bounds)),
            StageResult::Filtered { .. } => None,
        };

        Ok(PipelineResult {
            profile_name: prf.name.clone(),
            target_name: seq.name.clone(),
            profile_length: prf.length,
            target_length: seq.length,
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

pub fn run_pipeline_profile_to_sequence(pipeline: &mut Pipeline, seeds: Seeds2) {
    let tl_pipeline: ThreadLocal<RefCell<Pipeline>> = ThreadLocal::new();

    seeds
        .seeds
        .par_chunks(100)
        .panic_fuse()
        .try_for_each(|chunk| -> anyhow::Result<()> {
            let now = Instant::now();
            let mut pipeline = tl_pipeline
                .get_or(|| RefCell::new(pipeline.clone()))
                .borrow_mut();

            let res = chunk
                .iter()
                .map(|seed| pipeline.run(seed))
                .collect::<Result<Vec<_>, _>>()?;

            let output_stats = pipeline.output.run(&res)?;
            pipeline.stats.add_sample(&res, &output_stats);

            pipeline
                .stats
                .add_threaded_time(ThreadedTimed::Total, now.elapsed());

            Ok(())
        })
        .unwrap();
}
