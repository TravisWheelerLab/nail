use std::{
    io::{stdout, Write},
    time::{Duration, Instant},
};

use derive_builder::Builder;
use libnail::output::output_tabular::{Field, TableFormat};

use crate::{args::SearchArgs, pipeline::StageResult, util::PathExt};

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

#[derive(Builder, Default)]
#[builder(setter(strip_option), default)]
pub struct OutputStageStats {
    pub n_reported: usize,
    pub write_time: Duration,
}

impl OutputStageStatsBuilder {
    fn add_write_time(&mut self, duration: Duration) {
        match self.write_time {
            Some(ref mut time) => *time += duration,
            None => {
                self.write_time(duration);
            }
        }
    }
}

pub trait PipelineOutput: Send {
    fn output(&mut self, res: &PipelineResult) -> anyhow::Result<()>;
}

// ---

struct TableOutput<W>
where
    W: Write + Send,
{
    buf: W,
    table_format: TableFormat,
    header_written: bool,
}

impl<W> TableOutput<W>
where
    W: Write + Send,
{
    pub fn new(buf: W) -> Self {
        Self {
            buf,
            table_format: TableFormat::new(&DEFAULT_COLUMNS).expect("failed to build table format"),
            header_written: false,
        }
    }
}

impl<W> PipelineOutput for TableOutput<W>
where
    W: Write + Send,
{
    fn output(&mut self, res: &PipelineResult) -> anyhow::Result<()> {
        if !self.header_written {
            writeln!(self.buf, "{}", self.table_format.header()?)?;
            self.header_written = true;
        }

        // TODO: state for column width formatting
        if let Some(StageResult::Passed { data: ali, .. }) = &res.align_result {
            writeln!(self.buf, "{}", ali.tab_string_formatted(&self.table_format))?;
        }

        Ok(())
    }
}

// ---

struct AlignmentOutput<W>
where
    W: Write + Send,
{
    buf: W,
}

impl<W> AlignmentOutput<W>
where
    W: Write + Send,
{
    pub fn new(buf: W) -> Self {
        Self { buf }
    }
}

impl<W> PipelineOutput for AlignmentOutput<W>
where
    W: Write + Send,
{
    fn output(&mut self, res: &PipelineResult) -> anyhow::Result<()> {
        if let Some(StageResult::Passed { data: ali, .. }) = &res.align_result {
            writeln!(self.buf, "{}", ali.ali_string())?;
        }
        Ok(())
    }
}

// ---

pub struct OutputStage {
    output: Vec<Box<dyn PipelineOutput>>,
    e_value_threshold: f64,
}

impl OutputStage {
    pub fn new(args: &SearchArgs) -> anyhow::Result<Self> {
        let mut output: Vec<Box<dyn PipelineOutput>> = vec![];

        if args.ali_to_stdout {
            output.push(Box::new(AlignmentOutput::new(stdout())))
        } else if let Some(path) = &args.io_args.ali_results_path {
            output.push(Box::new(AlignmentOutput::new(path.open(true)?)))
        }

        if let Some(path) = &args.io_args.tbl_results_path {
            output.push(Box::new(TableOutput::new(path.open(true)?)));
        }

        Ok(Self {
            output,
            e_value_threshold: args.pipeline_args.e_value_threshold,
        })
    }

    pub fn run(&mut self, pipeline_results: &[PipelineResult]) -> anyhow::Result<OutputStageStats> {
        let mut stats = OutputStageStatsBuilder::default();

        let mut n_reported = 0;

        for res in pipeline_results {
            if let Some(StageResult::Passed { data: ali, .. }) = &res.align_result {
                if ali.scores.e_value <= self.e_value_threshold {
                    n_reported += 1;
                    for writer in &mut self.output {
                        let now = Instant::now();
                        writer.output(res)?;
                        stats.add_write_time(now.elapsed());
                    }
                }
            }
        }

        stats.n_reported(n_reported);

        Ok(stats.build()?)
    }
}
