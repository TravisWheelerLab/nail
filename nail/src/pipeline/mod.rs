mod cloud_stage;
pub use cloud_stage::*;

mod seed_stage;
pub use seed_stage::*;

mod align_stage;
pub use align_stage::*;

mod output_stage;
pub use output_stage::*;

use std::sync::Arc;
use std::sync::Mutex;

use libnail::{
    align::{structs::Seed, Bits},
    structs::Profile,
};

use crate::{io::Fasta, stats::Stats};

use anyhow::{bail, Context};

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
    /// The position of the profile among the queries
    pub profile_pos: usize,
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
    pub profiles: Arc<Vec<Profile>>,
    /// The profile of the last seed run, and its position.
    pub prf: Option<(usize, Profile)>,
    pub targets: Arc<Fasta>,
    /// Whether to check each seed's coordinates against its profile and target.
    pub check_seed_bounds: bool,
    pub cloud_search: Box<dyn CloudSearchStage>,
    pub align: Box<dyn AlignStage>,
    pub output: Arc<Mutex<OutputStage>>,
    pub stats: Stats,
}

/// Checks that a seed's 1-based, inclusive coordinates lie within
/// a profile of `prf_length` and a target of `seq_length`.
pub fn check_seed_bounds(seed: &Seed, prf_length: usize, seq_length: usize) -> anyhow::Result<()> {
    if seed.prf_start < 1 || seed.prf_end > prf_length {
        bail!(
            "query coordinates {}-{} are outside 1-{prf_length}",
            seed.prf_start,
            seed.prf_end
        );
    }
    if seed.prf_start > seed.prf_end {
        bail!(
            "query start {} is after query end {}",
            seed.prf_start,
            seed.prf_end
        );
    }
    if seed.seq_start < 1 || seed.seq_end > seq_length {
        bail!(
            "target coordinates {}-{} are outside 1-{seq_length}",
            seed.seq_start,
            seed.seq_end
        );
    }
    if seed.seq_start > seed.seq_end {
        bail!(
            "target start {} is after target end {}",
            seed.seq_start,
            seed.seq_end
        );
    }
    Ok(())
}

impl Pipeline {
    pub fn run(&mut self, seed: &Seed) -> anyhow::Result<PipelineResult> {
        let profile_is_cached = match &self.prf {
            Some((position, _)) => *position == seed.prf_pos,
            None => false,
        };

        if !profile_is_cached {
            let prf = self.profiles.get(seed.prf_pos).context("no prf")?;
            self.prf = Some((seed.prf_pos, prf.clone()));
        }

        let (_, prf) = self.prf.as_mut().expect("no prf");
        let seq = self.targets.get(seed.seq_pos).context("no seq")?;

        if self.check_seed_bounds {
            check_seed_bounds(seed, prf.length, seq.length).with_context(|| {
                format!(
                    "seed for query {} against target {} is out of range",
                    prf.name, seq.name
                )
            })?;
        }

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
            profile_pos: seed.prf_pos,
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

#[cfg(test)]
mod tests {
    use super::*;

    fn seed(prf_start: usize, prf_end: usize, seq_start: usize, seq_end: usize) -> Seed {
        Seed {
            prf_pos: 0,
            seq_pos: 0,
            prf_start,
            prf_end,
            seq_start,
            seq_end,
            score: 0.0,
            e_value: 0.0,
        }
    }

    #[test]
    fn a_seed_within_both_lengths_passes() -> anyhow::Result<()> {
        check_seed_bounds(&seed(1, 250, 1, 300), 250, 300)
    }

    #[test]
    fn a_seed_outside_either_length_fails() {
        let cases = [
            (
                seed(0, 10, 1, 10),
                "query coordinates 0-10 are outside 1-250",
            ),
            (
                seed(1, 251, 1, 10),
                "query coordinates 1-251 are outside 1-250",
            ),
            (seed(20, 10, 1, 10), "query start 20 is after query end 10"),
            (
                seed(1, 10, 0, 10),
                "target coordinates 0-10 are outside 1-300",
            ),
            (
                seed(1, 10, 1, 301),
                "target coordinates 1-301 are outside 1-300",
            ),
            (
                seed(1, 10, 20, 10),
                "target start 20 is after target end 10",
            ),
        ];

        for (seed, message) in cases {
            let Err(err) = check_seed_bounds(&seed, 250, 300) else {
                panic!("{seed:?} passed the bounds check");
            };
            assert_eq!(err.to_string(), message);
        }
    }
}
