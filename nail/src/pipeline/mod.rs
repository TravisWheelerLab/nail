mod cloud_step;
pub use cloud_step::*;

mod seed_step;
pub use seed_step::*;

mod align_step;
pub use align_step::*;

use std::collections::HashMap;
use std::io::{stdout, Write};
use std::sync::{Arc, Mutex};

use rayon::iter::{IntoParallelRefIterator, IntoParallelRefMutIterator, ParallelIterator};
use thiserror::Error;

use libnail::align::structs::Alignment;
use libnail::output::output_tabular::{Field, TableFormat};
use libnail::structs::{Hmm, Profile, Sequence};

use crate::args::AlignOutputArgs;
use crate::util::PathBufExt;

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

pub const DEFAULT_COLUMNS: [Field; 10] = [
    Field::Target,
    Field::Query,
    Field::TargetStart,
    Field::TargetEnd,
    Field::QueryStart,
    Field::QueryEnd,
    Field::Score,
    Field::CompBias,
    Field::Evalue,
    Field::CellFrac,
];

#[derive(Clone)]
pub struct Pipeline {
    pub seed: Box<dyn SeedStep + Send + Sync>,
    pub cloud_search: Box<dyn CloudSearchStep + Send + Sync>,
    pub align: Box<dyn AlignStep + Send + Sync>,
}

impl Pipeline {
    fn run(
        &mut self,
        profile: &mut Profile,
        targets: &HashMap<String, Sequence>,
    ) -> Option<Vec<Alignment>> {
        self.seed.run(profile, targets).map(|seeds| {
            seeds
                .iter()
                .filter_map(|(target_name, seed)| {
                    let target = targets.get(target_name).unwrap();
                    self.cloud_search
                        .run(profile, target, seed)
                        .and_then(|bounds| self.align.run(profile, target, bounds))
                })
                .collect()
        })
    }
}

pub enum HeaderStatus {
    Unwritten,
    Written,
}

pub struct Output {
    alignment_writer: Box<dyn Write + Send>,
    table_writer: Box<dyn Write + Send>,
    table_format: TableFormat,
    header_status: HeaderStatus,
}

impl Output {
    pub fn new(args: &AlignOutputArgs) -> anyhow::Result<Self> {
        Ok(Self {
            alignment_writer: match &args.ali_results_path {
                Some(path) => Box::new(path.open(true)?),
                None => Box::new(stdout()),
            },
            table_writer: Box::new(args.tsv_results_path.open(true)?),
            table_format: TableFormat::new(&DEFAULT_COLUMNS)?,
            header_status: HeaderStatus::Unwritten,
        })
    }

    pub fn write(&mut self, alignments: &mut [Alignment]) -> anyhow::Result<()> {
        self.table_format.reset_widths();
        self.table_format.update_widths(alignments);

        alignments.sort_by(|a, b| a.e_value.partial_cmp(&b.e_value).unwrap());

        if let HeaderStatus::Unwritten = self.header_status {
            let header = TableFormat::header(&self.table_format)?;
            writeln!(self.table_writer, "{header}")?;
            self.header_status = HeaderStatus::Written;
        }

        alignments.iter().for_each(|ali| {
            writeln!(
                self.table_writer,
                "{}",
                ali.tab_string_formatted(&self.table_format)
            )
            .expect("failed to write tabular output");

            writeln!(self.alignment_writer, "{}", ali.ali_string())
                .expect("failed to write alignment output");
        });

        Ok(())
    }
}

pub fn run_pipeline_profile_to_sequence(
    queries: &mut [Profile],
    targets: &HashMap<String, Sequence>,
    pipeline: Pipeline,
    output: Arc<Mutex<Output>>,
) {
    queries.par_iter_mut().for_each_with(
        (pipeline, targets, output),
        |(pipeline, targets, output), profile| {
            let mut alignments: Vec<Alignment> = pipeline.run(profile, targets).unwrap();

            match output.lock() {
                Ok(mut guard) => {
                    guard.write(&mut alignments).unwrap();
                }
                Err(_) => panic!(),
            }
        },
    )
}

pub fn run_pipeline_sequence_to_sequence(
    queries: &[Sequence],
    targets: &HashMap<String, Sequence>,
    pipeline: Pipeline,
    output: Arc<Mutex<Output>>,
) {
    queries.par_iter().for_each_with(
        (pipeline, targets, output),
        |(pipeline, targets, output), sequence| {
            let mut profile = Hmm::from_blosum_62_and_sequence(sequence)
                .map(|h| Profile::new(&h))
                .expect("failed to build profile from sequence");

            profile.calibrate_tau(200, 100, 0.04);

            let mut alignments: Vec<Alignment> = pipeline.run(&mut profile, targets).unwrap();

            match output.lock() {
                Ok(mut guard) => {
                    guard.write(&mut alignments).unwrap();
                }
                Err(_) => panic!(),
            }
        },
    )
}
