use std::{
    io::Write,
    time::{Duration, Instant},
};

use derive_builder::Builder;
use libnail::output::output_tabular::{Field, TableFormat};

use crate::pipeline::StageResult;

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

pub const BLAST_COLUMNS: [Field; 12] = [
    Field::Query,         //  1.  qseqid      query or source (gene) sequence id
    Field::Target,        //  2.  sseqid      subject or target (reference genome) sequence id
    Field::Pid,           //  3.  pident      percentage of identical positions
    Field::Length,        //  4.  length      alignment length (sequence overlap)
    Field::MismatchCount, //  5.  mismatch    number of mismatches
    Field::GapOpenCount,  //  6.  gapopen     number of gap openings
    Field::QueryStart,    //  7.  qstart      start of alignment in query
    Field::QueryEnd,      //  8.  qend        end of alignment in query
    Field::TargetStart,   //  9.  sstart      start of alignment in subject
    Field::TargetEnd,     // 10.  send        end of alignment in subject
    Field::Evalue,        // 11.  evalue      expect value
    Field::Score,         // 12.  bitscore    bit score
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

pub struct TableOutput<W>
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
    pub fn new(buf: W, format: TableFormat) -> Self {
        Self {
            buf,
            table_format: format,
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

pub struct AlignmentOutput<W>
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
    pub fn new(e_value_threshold: f64) -> anyhow::Result<Self> {
        Ok(Self {
            output: vec![],
            e_value_threshold,
        })
    }

    pub fn add<O>(&mut self, output: O)
    where
        O: PipelineOutput + 'static,
    {
        self.output.push(Box::new(output))
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
