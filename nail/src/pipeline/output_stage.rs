use std::{
    io::Write,
    sync::{Arc, Mutex},
    time::{Duration, Instant},
};

use anyhow::{anyhow, Context};
use derive_builder::Builder;
use libnail::{
    align::structs::Alignment,
    output::output_tabular::{Field, TableFormat},
};

use crate::{args::OutputArgs, util::PathBufExt};

use super::PipelineResult;

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
pub enum HeaderStatus {
    Unwritten,
    Written,
}

#[derive(Builder, Default)]
#[builder(setter(strip_option), default)]
pub struct OutputStageStats {
    pub lock_time: Duration,
    pub write_time: Duration,
}

impl OutputStageStatsBuilder {
    fn add_lock_time(&mut self, duration: Duration) {
        match self.lock_time {
            Some(ref mut time) => *time += duration,
            None => {
                self.lock_time(duration);
            }
        }
    }

    fn add_write_time(&mut self, duration: Duration) {
        match self.write_time {
            Some(ref mut time) => *time += duration,
            None => {
                self.write_time(duration);
            }
        }
    }
}

#[derive(Clone)]
pub struct OutputStage {
    alignment_writer: Option<Arc<Mutex<Box<dyn Write + Send>>>>,
    table_writer: Option<Arc<Mutex<Box<dyn Write + Send>>>>,
    stats_writer: Option<Arc<Mutex<Box<dyn Write + Send>>>>,
    e_value_threshold: f64,
    table_format: TableFormat,
    header_status: Arc<Mutex<HeaderStatus>>,
}

impl OutputStage {
    pub fn new(args: &OutputArgs) -> anyhow::Result<Self> {
        Ok(Self {
            alignment_writer: match &args.ali_results_path {
                Some(path) => Some(Arc::new(Mutex::new(Box::new(path.open(true)?)))),
                None => None,
            },
            table_writer: Some(Arc::new(Mutex::new(Box::new(
                args.tbl_results_path.open(true)?,
            )))),
            table_format: TableFormat::new(&DEFAULT_COLUMNS)?,
            e_value_threshold: args.e_value_threshold,
            header_status: Arc::new(Mutex::new(HeaderStatus::Unwritten)),
            stats_writer: match &args.stats_results_path {
                Some(path) => Some(Arc::new(Mutex::new(Box::new(path.open(true)?)))),
                None => None,
            },
        })
    }

    pub fn run(&mut self, pipeline_results: &[PipelineResult]) -> anyhow::Result<OutputStageStats> {
        let mut stats = OutputStageStatsBuilder::default();

        let mut reported: Vec<&Alignment> = pipeline_results
            .iter()
            .filter_map(|r| r.align_result.as_ref())
            .filter_map(|r| match r {
                super::StageResult::Filtered { .. } => None,
                super::StageResult::Passed { data, .. } => Some(data),
            })
            .filter(|a| a.scores.e_value <= self.e_value_threshold)
            .collect();

        reported.sort_by(|a, b| a.scores.e_value.partial_cmp(&b.scores.e_value).unwrap());

        if let Some(writer) = &self.alignment_writer {
            let now = Instant::now();
            match writer.lock() {
                Ok(mut guard) => {
                    stats.add_lock_time(now.elapsed());

                    let now = Instant::now();
                    reported
                        .iter()
                        .try_for_each(|ali| writeln!(guard, "{}", ali.ali_string()))
                        .with_context(|| "failed to write to alignment writer")?;

                    stats.add_write_time(now.elapsed());
                    Ok(())
                }
                Err(_) => Err(anyhow!("alignment writer mutex poisoned")),
            }?;
        }

        if let Some(writer) = &self.table_writer {
            let now = Instant::now();
            match writer.lock() {
                Ok(mut writer_guard) => {
                    stats.add_lock_time(now.elapsed());

                    self.table_format.reset_widths();
                    self.table_format.update_widths(&reported);

                    // TODO: it's a bit messy to put the header status in a
                    //       mutex I'd like to come up with a better way to
                    //       write the header just one time, but the problem
                    //       is that the output looks way better if we have
                    //       the first table format computed
                    match self.header_status.lock() {
                        Ok(mut header_status_guard) => {
                            if let HeaderStatus::Unwritten = *header_status_guard {
                                let header = TableFormat::header(&self.table_format)?;
                                writeln!(writer_guard, "{header}")?;
                                *header_status_guard = HeaderStatus::Written;
                            }
                            Ok(())
                        }
                        Err(_) => Err(anyhow!("header status mutex poisoned")),
                    }?;

                    let now = Instant::now();
                    reported.iter().try_for_each(|ali| {
                        writeln!(
                            writer_guard,
                            "{}",
                            ali.tab_string_formatted(&self.table_format)
                        )
                        .with_context(|| "failed to write to table writer")
                    })?;

                    stats.add_write_time(now.elapsed());
                    Ok(())
                }
                Err(_) => Err(anyhow!("table writer mutex poisoned")),
            }?;
        }

        if let Some(writer) = &self.stats_writer {
            let now = Instant::now();
            match writer.lock() {
                Ok(mut guard) => {
                    stats.add_lock_time(now.elapsed());

                    let now = Instant::now();
                    pipeline_results
                        .iter()
                        .try_for_each(|r| writeln!(guard, "{}", r.tab_string()))?;

                    stats.add_write_time(now.elapsed());
                    Ok(())
                }
                Err(_) => Err(anyhow!("stats writer mutex poisoned")),
            }?;
        }

        stats.build().map_err(Into::into)
    }
}
