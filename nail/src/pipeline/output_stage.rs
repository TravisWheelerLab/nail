use std::{
    io::Write,
    time::{Duration, Instant},
};

use derive_builder::Builder;
use libnail::align::{structs::Scores, Bits};
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

/// The test a hit must pass to be reported.
pub enum ReportThreshold {
    /// The largest E-value reported
    EValue(f64),
    /// The smallest bit score reported, one per profile, by profile position
    Score(Vec<Bits>),
}

impl ReportThreshold {
    /// Whether a hit of the profile at `profile_pos` with these scores is reported.
    pub fn passes(&self, profile_pos: usize, scores: &Scores) -> bool {
        match self {
            ReportThreshold::EValue(e_value) => scores.e_value <= *e_value,
            ReportThreshold::Score(cutoffs) => scores.bit_score >= cutoffs[profile_pos],
        }
    }
}

pub struct OutputStage {
    output: Vec<Box<dyn PipelineOutput>>,
    threshold: ReportThreshold,
}

impl OutputStage {
    pub fn new(threshold: ReportThreshold) -> anyhow::Result<Self> {
        Ok(Self {
            output: vec![],
            threshold,
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
                if self.threshold.passes(res.profile_pos, &ali.scores) {
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

#[cfg(test)]
mod tests {
    use super::*;

    fn scores(bit_score: f32, e_value: f64) -> Scores {
        Scores {
            forward_score: Bits(0.0),
            forward_p_value: 0.0,
            bit_score: Bits(bit_score),
            null_two_score: None,
            p_value: 0.0,
            e_value,
        }
    }

    #[test]
    fn score_threshold_compares_with_the_hits_own_profile() {
        let threshold = ReportThreshold::Score(vec![Bits(20.0), Bits(30.0)]);

        // the E-values are set to fail any E-value test,
        // so a pass here can only come from the score
        assert!(threshold.passes(1, &scores(30.0, 1e9)));
        assert!(!threshold.passes(1, &scores(29.9, 0.0)));
        assert!(threshold.passes(0, &scores(25.0, 1e9)));
    }

    #[test]
    fn e_value_threshold_reports_up_to_and_including_its_value() {
        let threshold = ReportThreshold::EValue(10.0);

        assert!(threshold.passes(0, &scores(0.0, 10.0)));
        assert!(!threshold.passes(0, &scores(100.0, 10.1)));
    }
}
