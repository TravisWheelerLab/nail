use libnail::structs::Sequence;
use std::{
    fmt::Debug,
    io::Write,
    sync::{
        atomic::{AtomicU64, Ordering},
        Arc,
    },
    time::Duration,
};
use strum::{EnumCount, EnumIter, IntoEnumIterator};

use crate::{
    pipeline::{
        OutputStageStats, PipelineResult,
        StageResult::{Filtered, Passed},
    },
    search::Queries,
};

#[repr(usize)]
#[derive(Clone, Copy, EnumIter, EnumCount)]
pub enum SerialTimed {
    Total,
    Seeding,
    Alignment,
}

#[repr(usize)]
#[derive(Clone, Copy, EnumIter, EnumCount)]
pub enum ThreadedTimed {
    Total,
    MemoryInit,
    HmmBuild,
    CloudSearch,
    Forward,
    Backward,
    Posterior,
    Traceback,
    NullTwo,
    OutputWrite,
    OutputMutex,
}

impl Debug for ThreadedTimed {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let str = match self {
            ThreadedTimed::Total => "total",
            ThreadedTimed::OutputWrite => "output (write)",
            ThreadedTimed::OutputMutex => "output (mutex)",
            ThreadedTimed::MemoryInit => "memory init",
            ThreadedTimed::HmmBuild => "hmm build",
            ThreadedTimed::CloudSearch => "cloud search",
            ThreadedTimed::Forward => "forward",
            ThreadedTimed::Backward => "backward",
            ThreadedTimed::Posterior => "posterior",
            ThreadedTimed::Traceback => "traceback",
            ThreadedTimed::NullTwo => "null two",
        };

        write!(f, "{}", str)
    }
}

#[repr(usize)]
#[derive(Clone, Copy, EnumIter, EnumCount)]
pub enum ComputedValue {
    Queries,
    Targets,
    Alignments,
    Cells,
}

impl Debug for ComputedValue {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let str = match self {
            ComputedValue::Queries => "queries",
            ComputedValue::Targets => "targets",
            ComputedValue::Alignments => "total alignments",
            ComputedValue::Cells => "total cells",
        };

        write!(f, "{}", str)
    }
}

#[repr(usize)]
#[derive(Clone, Copy, EnumIter, EnumCount)]
pub enum CountedValue {
    Seeds,
    PassedCloud,
    PassedForward,
    SeedCells,
    CloudForwardCells,
    CloudBackwardCells,
    ForwardCells,
    BackwardCells,
}

impl Debug for CountedValue {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let str = match self {
            CountedValue::Seeds => "seeds",
            CountedValue::PassedCloud => "passed cloud filter",
            CountedValue::PassedForward => "passed forward filter",
            CountedValue::SeedCells => "seed cells",
            CountedValue::CloudForwardCells => "cloud forward cells computed",
            CountedValue::CloudBackwardCells => "cloud backward cells computed",
            CountedValue::ForwardCells => "forward cells computed",
            CountedValue::BackwardCells => "backward cells computed",
        };

        write!(f, "{}", str)
    }
}

#[derive(Clone, Default)]
pub struct Stats {
    serial_times: [Duration; SerialTimed::COUNT],
    threaded_times: Arc<[AtomicU64; ThreadedTimed::COUNT]>,
    threaded_times_num_samples: Arc<[AtomicU64; ThreadedTimed::COUNT]>,
    counted_values: Arc<[AtomicU64; CountedValue::COUNT]>,
    computed_values: [u64; ComputedValue::COUNT],
}

impl Stats {
    pub fn new(queries: &Queries, targets: &[Sequence]) -> Self {
        let mut stats = Self::default();

        stats.set_computed_value(ComputedValue::Queries, queries.len() as u64);
        stats.set_computed_value(ComputedValue::Targets, targets.len() as u64);
        stats.set_computed_value(
            ComputedValue::Alignments,
            (queries.len() * targets.len()) as u64,
        );
        stats.set_computed_value(
            ComputedValue::Cells,
            queries
                .lengths()
                .iter()
                .flat_map(|len| targets.iter().map(move |t| (t.length * len) as u64))
                .sum(),
        );

        stats
    }

    pub fn add_sample(
        &mut self,
        pipeline_results: &[PipelineResult],
        output_stats: &OutputStageStats,
    ) {
        pipeline_results.iter().for_each(|result| {
            self.increment_count(CountedValue::Seeds);
            self.add_count(
                CountedValue::SeedCells,
                result.profile_length * result.target_length,
            );

            if let Some(ref result) = result.cloud_result {
                match result {
                    Filtered { stats } => {
                        self.add_count(CountedValue::CloudForwardCells, stats.forward_cells);
                        self.add_count(CountedValue::CloudBackwardCells, stats.backward_cells);

                        self.add_threaded_time(ThreadedTimed::CloudSearch, stats.memory_init_time);
                        self.add_threaded_time(ThreadedTimed::CloudSearch, stats.forward_time);
                        self.add_threaded_time(ThreadedTimed::CloudSearch, stats.backward_time);
                    }
                    Passed { stats, .. } => {
                        self.increment_count(CountedValue::PassedCloud);

                        self.add_count(CountedValue::CloudForwardCells, stats.forward_cells);
                        self.add_count(CountedValue::CloudBackwardCells, stats.backward_cells);

                        self.add_threaded_time(ThreadedTimed::CloudSearch, stats.memory_init_time);
                        self.add_threaded_time(ThreadedTimed::CloudSearch, stats.forward_time);
                        self.add_threaded_time(ThreadedTimed::CloudSearch, stats.backward_time);

                        self.add_threaded_time(ThreadedTimed::CloudSearch, stats.reorient_time);
                        self.add_threaded_time(ThreadedTimed::CloudSearch, stats.merge_time);
                        self.add_threaded_time(ThreadedTimed::CloudSearch, stats.trim_time);
                    }
                }
            }

            if let Some(ref result) = result.align_result {
                match result {
                    Filtered { stats } => {
                        self.add_count(CountedValue::ForwardCells, stats.forward_cells);

                        self.add_threaded_time(ThreadedTimed::MemoryInit, stats.memory_init_time);
                        self.add_threaded_time(ThreadedTimed::Forward, stats.forward_time);
                    }
                    Passed { stats, .. } => {
                        self.increment_count(CountedValue::PassedForward);

                        self.add_count(CountedValue::ForwardCells, stats.forward_cells);
                        self.add_count(CountedValue::BackwardCells, stats.backward_cells);

                        self.add_threaded_time(ThreadedTimed::MemoryInit, stats.memory_init_time);
                        self.add_threaded_time(ThreadedTimed::Forward, stats.forward_time);

                        self.add_threaded_time(ThreadedTimed::Backward, stats.backward_time);
                        self.add_threaded_time(ThreadedTimed::Posterior, stats.posterior_time);
                        self.add_threaded_time(ThreadedTimed::Traceback, stats.traceback_time);
                        self.add_threaded_time(ThreadedTimed::NullTwo, stats.null_two_time);
                    }
                }
            }
        });
    }

    pub fn set_serial_time(&mut self, timed: SerialTimed, time: Duration) {
        self.serial_times[timed as usize] = time;
    }

    pub fn add_threaded_time(&mut self, timed: ThreadedTimed, time: Duration) {
        let time_nanos = Self::nanos(time);
        self.threaded_times[timed as usize].fetch_add(time_nanos, Ordering::SeqCst);
        self.threaded_times_num_samples[timed as usize].fetch_add(1, Ordering::SeqCst);
    }

    fn serial_time_total(&self, timed: SerialTimed) -> Duration {
        self.serial_times[timed as usize]
    }

    fn threaded_time_total(&self, timed: ThreadedTimed) -> Duration {
        let nanos = self.threaded_times[timed as usize].load(Ordering::SeqCst);
        Duration::from_nanos(nanos)
    }

    fn serial_time_pct(&self, timed: SerialTimed) -> f64 {
        let total_nanos = Self::nanos(self.serial_times[ThreadedTimed::Total as usize]) as f64;
        let nanos = Self::nanos(self.serial_times[timed as usize]) as f64;

        nanos / total_nanos
    }

    fn threaded_time_pct(&self, timed: ThreadedTimed) -> f64 {
        let total_nanos = Self::nanos(self.threaded_time_total(ThreadedTimed::Total)) as f64;
        let nanos = Self::nanos(self.threaded_time_total(timed)) as f64;

        nanos / total_nanos
    }

    fn threaded_untimed_pct(&self) -> f64 {
        let total_nanos = Self::nanos(self.threaded_time_total(ThreadedTimed::Total)) as f64;
        let untimed_nanos = Self::nanos(self.threaded_untimed_total()) as f64;

        untimed_nanos / total_nanos
    }

    fn threaded_untimed_total(&self) -> Duration {
        let total = self.threaded_time_total(ThreadedTimed::Total);

        let timed_sum = self.threaded_times[1..]
            .iter()
            .map(|t| Duration::from_nanos(t.load(Ordering::SeqCst)))
            .sum();

        total - timed_sum
    }

    fn computed_value(&self, computed: ComputedValue) -> u64 {
        self.computed_values[computed as usize]
    }

    fn set_computed_value(&mut self, computed: ComputedValue, value: u64) {
        self.computed_values[computed as usize] = value
    }

    fn counted_value(&self, counted: CountedValue) -> u64 {
        self.counted_values[counted as usize].load(Ordering::SeqCst)
    }

    pub fn increment_count(&mut self, counted: CountedValue) {
        self.counted_values[counted as usize].fetch_add(1, Ordering::SeqCst);
    }

    pub fn add_count(&mut self, counted: CountedValue, count: usize) {
        self.counted_values[counted as usize].fetch_add(count as u64, Ordering::SeqCst);
    }

    pub fn serial_string(&self, timed: SerialTimed) -> String {
        let width = format!(
            "{:.2}",
            self.serial_time_total(SerialTimed::Total).as_secs_f64()
        )
        .len();

        format!(
            "{:w$.2}s ({:>5.2}%)",
            self.serial_time_total(timed).as_secs_f64(),
            self.serial_time_pct(timed) * 100.0,
            w = width,
        )
    }

    pub fn write(&self, out: &mut impl Write) -> anyhow::Result<()> {
        self.write_stats(out)?;
        writeln!(out)?;
        self.write_runtime(out)
    }

    pub fn write_stats(&self, out: &mut impl Write) -> anyhow::Result<()> {
        let max_width = ComputedValue::iter()
            .map(|c| format!("{c:?}: {}", Self::format_num(self.computed_value(c))).len())
            .chain(
                CountedValue::iter()
                    .map(|c| format!("{c:?}: {}", Self::format_num(self.counted_value(c))).len()),
            )
            .max()
            .unwrap_or(0);

        ComputedValue::iter().try_for_each(|c| {
            let label = format!("{c:?}");
            let label_width = label.len();
            let count = Self::format_num(self.computed_value(c));
            writeln!(out, "{label}: {count:>w$}", w = max_width - label_width)
        })?;

        CountedValue::iter().try_for_each(|c| {
            let label = format!("{c:?}");
            let label_width = label.len();
            let count = Self::format_num(self.counted_value(c));
            writeln!(out, "{label}: {count:>w$}", w = max_width - label_width)
        })?;

        Ok(())
    }

    pub fn write_runtime(&self, out: &mut impl Write) -> anyhow::Result<()> {
        writeln!(out, "runtime: {}", self.serial_string(SerialTimed::Total),)?;

        writeln!(
            out,
            " ├─ seeding:   {}",
            self.serial_string(SerialTimed::Seeding)
        )?;

        writeln!(
            out,
            " └─ alignment: {}",
            self.serial_string(SerialTimed::Alignment)
        )?;

        let max_width = ThreadedTimed::iter()
            .map(|t| format!("{t:?}: {:.2}", self.threaded_time_total(t).as_secs_f64()).len())
            .max()
            .unwrap_or(0);

        ThreadedTimed::iter()
            .skip(1)
            .filter(|t| !self.threaded_time_total(*t).is_zero())
            .try_for_each(|timed| {
                let label_width = format!("{timed:?}").len();

                writeln!(
                    out,
                    "     ├─ {timed:?}: {:>w$.2}s ({:5.2}%)",
                    self.threaded_time_total(timed).as_secs_f64(),
                    self.threaded_time_pct(timed) * 100.0,
                    w = max_width - label_width
                )
            })?;

        let last_label = "[???]";
        let last_label_width = last_label.len();

        writeln!(
            out,
            "     └─ {last_label}: {:>w$.2}s ({:5.2}%)",
            self.threaded_untimed_total().as_secs_f64(),
            self.threaded_untimed_pct() * 100.0,
            w = max_width - last_label_width
        )?;

        Ok(())
    }

    pub fn nanos(time: Duration) -> u64 {
        // u64::MAX nanoseconds is like 5,000,000 hours
        // or something, so this clamp should be fine.
        time.as_nanos().min(u64::MAX as u128) as u64
    }

    pub fn format_num(num: u64) -> String {
        let num_str = num.to_string();
        let mut result = String::new();
        let len = num_str.len();

        for (i, ch) in num_str.chars().enumerate() {
            if i > 0 && (len - i) % 3 == 0 {
                result.push(',');
            }
            result.push(ch);
        }
        result
    }
}
