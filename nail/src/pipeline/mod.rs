mod cloud_step;
pub use cloud_step::*;

mod seed_step;
pub use seed_step::*;

mod align_step;
pub use align_step::*;

mod output_step;
pub use output_step::*;

use std::collections::HashMap;
use std::io::stdout;
use std::sync::{Arc, Mutex};
use std::time::Instant;

use rayon::iter::{IntoParallelRefIterator, IntoParallelRefMutIterator, ParallelIterator};
use thiserror::Error;

use libnail::align::structs::{Alignment, RowBounds, Seed};
use libnail::structs::{Hmm, Profile, Sequence};

use crate::stats::Stats;

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

#[derive(Clone)]
pub struct Pipeline {
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
    ) -> Option<Vec<Alignment>> {
        let seeds = self.seed.run(profile, targets)?;

        let mut alignments: Vec<Alignment> = seeds
            .iter()
            .filter_map(|(target_name, seed)| {
                let target = targets.get(target_name)?;
                self.stats.increment_samples();

                let now = Instant::now();
                let maybe_bounds = self.cloud_search.run(profile, target, seed);
                self.stats.add_cloud_search_time(now.elapsed());

                let bounds = match maybe_bounds {
                    Some(bounds) => {
                        self.stats.increment_passed_cloud();
                        bounds
                    }
                    None => return None,
                };

                let now = Instant::now();
                let maybe_ali = self.align.run(profile, target, bounds);
                self.stats.add_forward_time(now.elapsed());

                match maybe_ali {
                    Some(ali) => {
                        self.stats.increment_passed_forward();
                        Some(ali)
                    }
                    None => None,
                }
            })
            .collect();

        let now = Instant::now();
        match self.output.lock() {
            Ok(mut guard) => guard.write(&mut alignments).unwrap(),
            Err(_) => panic!(),
        };
        self.stats.add_output_time(now.elapsed());

        Some(alignments)
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
            pipeline.run(profile, targets);
        },
    );

    pipeline.stats.write(&mut stdout());
}

pub fn run_pipeline_sequence_to_sequence(
    queries: &[Sequence],
    targets: &HashMap<String, Sequence>,
    pipeline: Pipeline,
) {
    queries.par_iter().panic_fuse().for_each_with(
        (pipeline, targets),
        |(pipeline, targets), sequence| {
            let mut profile = Hmm::from_blosum_62_and_sequence(sequence)
                .map(|h| Profile::new(&h))
                .expect("failed to build profile from sequence");

            profile.calibrate_tau(200, 100, 0.04);

            pipeline.run(&mut profile, targets).unwrap();
        },
    )
}
