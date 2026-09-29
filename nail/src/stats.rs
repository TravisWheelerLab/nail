use std::{
    collections::HashMap,
    fmt::Debug,
    io::Write,
    sync::{
        atomic::{AtomicU64, AtomicUsize, Ordering},
        Arc,
    },
    time::Duration,
};

use anyhow::{anyhow, Context};
use strum::{EnumCount, EnumIter, IntoEnumIterator};

use crate::{
    args::SearchArgs,
    pipeline::{
        OutputStageStats, PipelineResult,
        StageResult::{Filtered, Passed},
    },
    search::Queries,
    util::PathExt,
};

pub struct Bytes(pub usize);

#[allow(dead_code)]
impl Bytes {
    pub fn kib(&self) -> f64 {
        self.0 as f64 / 2.0_f64.powi(10)
    }

    pub fn mib(&self) -> f64 {
        self.0 as f64 / 2.0_f64.powi(20)
    }

    pub fn gib(&self) -> f64 {
        self.0 as f64 / 2.0_f64.powi(30)
    }
}

impl std::ops::Add for Bytes {
    type Output = Bytes;

    fn add(self, rhs: Self) -> Self::Output {
        Bytes(self.0 + rhs.0)
    }
}

impl std::ops::Sub for Bytes {
    type Output = Bytes;

    fn sub(self, rhs: Self) -> Self::Output {
        Bytes(self.0 - rhs.0)
    }
}

impl std::iter::Sum for Bytes {
    fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
        iter.fold(Bytes(0), std::ops::Add::add)
    }
}

#[repr(usize)]
#[derive(Clone, Copy, EnumIter, EnumCount)]
pub enum SerialTimed {
    Total,
    Setup,
    Seeding,
    Alignment,
}

impl Debug for SerialTimed {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let str = match self {
            SerialTimed::Total => "total",
            SerialTimed::Setup => "setup",
            SerialTimed::Seeding => "seeding",
            SerialTimed::Alignment => "alignment",
        };

        write!(f, "{}", str)
    }
}

#[repr(usize)]
#[derive(Clone, Copy, EnumIter, EnumCount)]
pub enum SetupTimed {
    Total,
    QueryIndex,
    TargetIndex,
    PipelineBuild,
}

impl Debug for SetupTimed {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let str = match self {
            SetupTimed::Total => "total",
            SetupTimed::QueryIndex => "query index",
            SetupTimed::TargetIndex => "target index",
            SetupTimed::PipelineBuild => "pipeline build",
        };

        write!(f, "{}", str)
    }
}

#[repr(usize)]
#[derive(Clone, Copy, EnumIter, EnumCount)]
pub enum SeedTimed {
    Total,
    DbWrite,
    Prefilter,
    Slice,
    Align,
    Decide,
    Merge,
    Convertalis,
    Index,
}

impl Debug for SeedTimed {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let str = match self {
            SeedTimed::Total => "total",
            SeedTimed::DbWrite => "db write",
            SeedTimed::Prefilter => "prefilter",
            SeedTimed::Slice => "slice",
            SeedTimed::Align => "align",
            SeedTimed::Decide => "decide",
            SeedTimed::Merge => "merge",
            SeedTimed::Convertalis => "convertalis",
            SeedTimed::Index => "seed index",
        };

        write!(f, "{}", str)
    }
}

#[repr(usize)]
#[derive(Clone, Copy, EnumIter, EnumCount)]
pub enum ThreadedTimed {
    Total,
    MemoryInit,
    CloudSearch,
    Forward,
    Backward,
    Posterior,
    OptimalAccuracy,
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
            ThreadedTimed::CloudSearch => "cloud search",
            ThreadedTimed::Forward => "forward",
            ThreadedTimed::Backward => "backward",
            ThreadedTimed::Posterior => "posterior",
            ThreadedTimed::OptimalAccuracy => "optimal accuracy",
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
            ComputedValue::Alignments => "total potential alignments",
            ComputedValue::Cells => "total potential DP cells",
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
    PassedReport,
    SeedCells,
    CloudForwardCells,
    CloudBackwardCells,
    ForwardCells,
    BackwardCells,
}

impl Debug for CountedValue {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let str = match self {
            CountedValue::Seeds => "passed seed filter",
            CountedValue::PassedCloud => "passed cloud filter",
            CountedValue::PassedForward => "passed forward filter",
            CountedValue::PassedReport => "passed reporting filter",
            CountedValue::SeedCells => "total potential seed DP cells",
            CountedValue::CloudForwardCells => "cloud forward DP cells computed",
            CountedValue::CloudBackwardCells => "cloud backward DP cells computed",
            CountedValue::ForwardCells => "forward DP cells computed",
            CountedValue::BackwardCells => "backward DP cells computed",
        };

        write!(f, "{}", str)
    }
}

#[derive(Clone, Default)]
pub struct Stats {
    serial_times: [Duration; SerialTimed::COUNT],
    setup_times: [Duration; SetupTimed::COUNT],
    seed_times: [Duration; SeedTimed::COUNT],
    threaded_times: Arc<[AtomicU64; ThreadedTimed::COUNT]>,
    threaded_times_num_samples: Arc<[AtomicU64; ThreadedTimed::COUNT]>,
    counted_values: Arc<[AtomicU64; CountedValue::COUNT]>,
    computed_values: [u64; ComputedValue::COUNT],
    seed_counts_by_query: HashMap<String, u64>,
    hit_counts_by_query: Arc<HashMap<String, AtomicUsize>>,
    num_threads: usize,
}

impl Stats {
    pub fn new(queries: &Queries, n_targets: usize) -> Self {
        let mut stats = Self::default();

        let (query_names, n_queries) = match queries {
            Queries::Sequence(db) => (db.index.keys(), db.len()),
            Queries::Profile(db) => (db.index.keys(), db.len()),
        };

        stats.hit_counts_by_query = Arc::new(
            query_names
                .map(|n| (n.to_string(), AtomicUsize::new(0)))
                .collect::<HashMap<_, _>>(),
        );

        stats.set_computed_value(ComputedValue::Queries, n_queries as u64);
        stats.set_computed_value(ComputedValue::Targets, n_targets as u64);
        stats.set_computed_value(ComputedValue::Alignments, (n_queries * n_targets) as u64);

        stats
    }

    pub fn set_num_threads(&mut self, num_threads: usize) {
        self.num_threads = num_threads;
    }

    pub fn set_setup_time(&mut self, timed: SetupTimed, time: Duration) {
        self.setup_times[timed as usize] = time;
    }

    pub fn set_seed_time(&mut self, timed: SeedTimed, time: Duration) {
        self.seed_times[timed as usize] = time;
    }

    pub fn add_seed_time(&mut self, timed: SeedTimed, time: Duration) {
        self.seed_times[timed as usize] += time;
    }

    pub fn set_seed_counts(&mut self, counts: HashMap<String, u64>) {
        self.seed_counts_by_query = counts;
    }

    pub fn add_sample(
        &mut self,
        pipeline_results: &[PipelineResult],
        output_stats: &OutputStageStats,
    ) {
        pipeline_results.iter().for_each(|pl_result| {
            self.increment_count(CountedValue::Seeds);
            self.add_count(
                CountedValue::SeedCells,
                pl_result.profile_length * pl_result.target_length,
            );

            if let Some(ref cld_result) = pl_result.cloud_result {
                match cld_result {
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

            if let Some(ref al_result) = pl_result.align_result {
                match al_result {
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
                        self.add_threaded_time(
                            ThreadedTimed::OptimalAccuracy,
                            stats.optimal_accuracy_time,
                        );
                        self.add_threaded_time(ThreadedTimed::Traceback, stats.traceback_time);
                        self.add_threaded_time(ThreadedTimed::NullTwo, stats.null_two_time);

                        self.hit_counts_by_query[&pl_result.profile_name]
                            .fetch_add(1, Ordering::Relaxed);
                    }
                }
            }
        });
        self.add_count(CountedValue::PassedReport, output_stats.n_reported);
        self.add_threaded_time(ThreadedTimed::OutputWrite, output_stats.write_time);
    }

    pub fn set_serial_time(&mut self, timed: SerialTimed, time: Duration) {
        self.serial_times[timed as usize] = time;
    }

    pub fn add_threaded_time(&mut self, timed: ThreadedTimed, time: Duration) {
        let time_nanos = Self::nanos(time);
        self.threaded_times[timed as usize].fetch_add(time_nanos, Ordering::SeqCst);
        self.threaded_times_num_samples[timed as usize].fetch_add(1, Ordering::SeqCst);
    }

    fn seed_time_total(&self, timed: SeedTimed) -> Duration {
        self.seed_times[timed as usize]
    }

    fn serial_time_total(&self, timed: SerialTimed) -> Duration {
        self.serial_times[timed as usize]
    }

    fn threaded_time_total(&self, timed: ThreadedTimed) -> Duration {
        let nanos = self.threaded_times[timed as usize].load(Ordering::SeqCst);
        Duration::from_nanos(nanos)
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

    fn serial_untimed_total(&self) -> Duration {
        let total = self.serial_time_total(SerialTimed::Total);
        let timed_sum = self.serial_times[1..].iter().sum();

        total.saturating_sub(timed_sum)
    }

    pub fn serial_string(&self, timed: SerialTimed) -> String {
        self.serial_duration_string(self.serial_time_total(timed))
    }

    fn serial_duration_string(&self, time: Duration) -> String {
        let total = self.serial_time_total(SerialTimed::Total);
        let width = format!("{:.2}", total.as_secs_f64()).len();

        format!(
            "{:w$.2}s ({:>5.2}%)",
            time.as_secs_f64(),
            Self::pct(time, total) * 100.0,
            w = width,
        )
    }

    fn pct(part: Duration, total: Duration) -> f64 {
        if total.is_zero() {
            0.0
        } else {
            part.as_secs_f64() / total.as_secs_f64()
        }
    }

    pub fn write_max_seqs_report(&self, args: &SearchArgs) -> anyhow::Result<()> {
        let queries: Vec<&String> = self.seed_counts_by_query.keys().collect();

        let path = args
            .io_args
            .tmp_dir_path
            .as_ref()
            .context("args.io_args.tmp_dir_path is somehow unset")?
            .join("max-seqs-report.txt");

        let mut out = path.open(true)?;

        let mut recs: Vec<(&String, u64, u64)> = queries
            .iter()
            .map(|&q| {
                let a = *self.seed_counts_by_query.get(q).unwrap_or(&0);
                let b = self
                    .hit_counts_by_query
                    .get(q)
                    .map_or(0u64, |a| a.load(Ordering::Relaxed) as u64);
                (q, a, b)
            })
            .collect();

        recs.sort_by(|a, b| b.1.cmp(&a.1));

        use crate::util::term::*;
        if let Some(first) = recs.first() {
            if first.1 == args.seed_args.max_seqs as u64 {
                println!();
                println!(
                "{RED}warning{RESET}: one or more queries saturated the mmseqs {YELLOW}--max-seqs{RESET} limit of {YELLOW}{}{RESET}",
                args.seed_args.max_seqs
            );
                println!("         for a full report, view the file: {YELLOW}{path:?}{RESET}",);
            }
        }

        let h1 = "query";
        let h2 = "seeds";
        let h3 = "reported hits";

        let w1 = h1
            .len()
            .max(queries.iter().map(|q| q.len()).max().unwrap_or(0));
        let w2 = h2.len().max(args.seed_args.max_seqs.to_string().len());
        let w3 = h3.len().max(args.seed_args.max_seqs.to_string().len());

        writeln!(out, "{h1:W1$} {h2:W2$} {h3:W3$}", W1 = w1, W2 = w2, W3 = w3)?;
        writeln!(
            out,
            "{:W1$} {:W2$} {:W3$}",
            "-".repeat(w1),
            "-".repeat(w2),
            "-".repeat(w3),
            W1 = w1,
            W2 = w2,
            W3 = w3
        )?;

        recs.into_iter().try_for_each(|(q, a, b)| {
            writeln!(out, "{q:W1$} {a:W2$} {b:W3$}", W1 = w1, W2 = w2, W3 = w3)
        })?;

        Ok(())
    }

    pub fn write(&self, out: &mut impl Write) -> anyhow::Result<()> {
        writeln!(out)?;
        writeln!(out, "summary statistics:")?;
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
            if self.computed_value(c) > 0 {
                let label = format!("{c:?}");
                let label_width = label.len();
                let count = Self::format_num(self.computed_value(c));
                writeln!(out, " ├─ {label}: {count:>w$}", w = max_width - label_width)
            } else {
                Ok(())
            }
        })?;

        let values: Vec<_> = CountedValue::iter().collect();
        values.iter().take(values.len() - 1).try_for_each(|c| {
            let label = format!("{c:?}");
            let label_width = label.len();
            let count = Self::format_num(self.counted_value(*c));
            writeln!(out, " ├─ {label}: {count:>w$}", w = max_width - label_width)
        })?;

        let last = values
            .last()
            .ok_or(anyhow!("no CountedValues in Stats::write_stats()"))?;
        let label = format!("{last:?}");
        let label_width = label.len();
        let count = Self::format_num(self.counted_value(*last));
        writeln!(out, " └─ {label}: {count:>w$}", w = max_width - label_width)?;

        Ok(())
    }

    pub fn write_runtime(&self, out: &mut impl Write) -> anyhow::Result<()> {
        writeln!(out, "runtime: {}", self.serial_string(SerialTimed::Total))?;

        let misc = "[misc.]";
        let branch_width = SerialTimed::iter()
            .map(|t| format!("{t:?}").len())
            .max()
            .unwrap_or(0)
            .max(misc.len())
            + 1;

        let branch = |out: &mut dyn Write, label: String, time: String| {
            writeln!(out, " └─ {label:<branch_width$} {time}")
        };

        branch(
            out,
            format!("{:?}:", SerialTimed::Setup),
            self.serial_string(SerialTimed::Setup),
        )?;
        Self::write_leaves(
            out,
            SetupTimed::iter()
                .skip(1)
                .map(|t| (format!("{t:?}"), self.setup_times[t as usize])),
            self.setup_times[SetupTimed::Total as usize],
        )?;

        branch(
            out,
            format!("{:?}:", SerialTimed::Seeding),
            self.serial_string(SerialTimed::Seeding),
        )?;
        Self::write_leaves(
            out,
            SeedTimed::iter()
                .skip(1)
                .map(|t| (format!("{t:?}"), self.seed_time_total(t))),
            self.seed_time_total(SeedTimed::Total),
        )?;

        let wall = self.serial_time_total(SerialTimed::Alignment);
        let cpu = self.threaded_time_total(ThreadedTimed::Total);
        let busy = Self::pct(cpu, wall * self.num_threads as u32) * 100.0;
        branch(
            out,
            format!("{:?}:", SerialTimed::Alignment),
            format!(
                "{}   [{} threads, cpu {}, {busy:.1}% busy]",
                self.serial_string(SerialTimed::Alignment),
                self.num_threads,
                Self::format_secs(cpu),
            ),
        )?;
        self.write_alignment_leaves(out, wall, cpu)?;

        branch(
            out,
            format!("{misc}:"),
            self.serial_duration_string(self.serial_untimed_total()),
        )?;

        Ok(())
    }

    fn write_alignment_leaves(
        &self,
        out: &mut impl Write,
        wall: Duration,
        cpu: Duration,
    ) -> anyhow::Result<()> {
        let mut rows: Vec<(String, Duration)> = ThreadedTimed::iter()
            .skip(1)
            .map(|t| (format!("{t:?}"), self.threaded_time_total(t)))
            .filter(|(_, t)| !t.is_zero())
            .collect();
        let timed_sum: Duration = rows.iter().map(|(_, t)| *t).sum();
        rows.push(("[misc.]".to_string(), cpu.saturating_sub(timed_sum)));

        let rows: Vec<(String, String, String, String)> = rows
            .into_iter()
            .map(|(label, leaf_cpu)| {
                let share = Self::pct(leaf_cpu, cpu);
                (
                    label,
                    format!("{:.2}s", wall.as_secs_f64() * share),
                    format!("{:.2}%", share * 100.0),
                    Self::format_secs(leaf_cpu),
                )
            })
            .collect();

        let width = |header: &str, col: fn(&(String, String, String, String)) -> &String| {
            rows.iter()
                .map(|r| col(r).len())
                .max()
                .unwrap_or(0)
                .max(header.len())
        };
        let label_w = width("", |r| &r.0);
        let wall_w = width("wall", |r| &r.1);
        let pct_w = width("%", |r| &r.2);
        let cpu_w = width("cpu", |r| &r.3);

        writeln!(
            out,
            "     │  {:label_w$}   {:>wall_w$}   {:>pct_w$}   {:>cpu_w$}",
            "", "wall", "%", "cpu"
        )?;
        writeln!(
            out,
            "     │  {:label_w$}   {}   {}   {}",
            "",
            "-".repeat(wall_w),
            "-".repeat(pct_w),
            "-".repeat(cpu_w)
        )?;

        let last = rows.len() - 1;
        rows.iter()
            .enumerate()
            .try_for_each(|(i, (label, wall, pct, cpu))| {
                let glyph = if i == last { "└─" } else { "├─" };
                writeln!(
                    out,
                    "     {glyph} {label:label_w$}   {wall:>wall_w$}   {pct:>pct_w$}   {cpu:>cpu_w$}"
                )
            })?;

        Ok(())
    }

    fn format_secs(time: Duration) -> String {
        let secs = time.as_secs_f64();
        let whole = Self::format_num(secs.trunc() as u64);
        let frac = (secs.fract() * 100.0).round() as u64;
        format!("{whole}.{frac:02}s")
    }

    fn write_leaves(
        out: &mut impl Write,
        leaves: impl Iterator<Item = (String, Duration)>,
        total: Duration,
    ) -> anyhow::Result<()> {
        let mut rows: Vec<(String, Duration)> = leaves.filter(|(_, t)| !t.is_zero()).collect();
        let timed_sum: Duration = rows.iter().map(|(_, t)| *t).sum();
        rows.push(("[misc.]".to_string(), total.saturating_sub(timed_sum)));

        let max_width = rows
            .iter()
            .map(|(label, t)| format!("{label}: {:.2}", t.as_secs_f64()).len())
            .max()
            .unwrap_or(0);

        let last = rows.len() - 1;
        rows.iter().enumerate().try_for_each(|(i, (label, t))| {
            let glyph = if i == last { "└─" } else { "├─" };
            writeln!(
                out,
                "     {glyph} {label}: {:>w$.2}s ({:5.2}%)",
                t.as_secs_f64(),
                Self::pct(*t, total) * 100.0,
                w = max_width - label.len()
            )
        })?;

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
            if i > 0 && (len - i).is_multiple_of(3) {
                result.push(',');
            }
            result.push(ch);
        }
        result
    }
}
