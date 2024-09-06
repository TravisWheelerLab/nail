use std::{
    fmt::Display,
    io::Write,
    sync::{
        atomic::{AtomicUsize, Ordering},
        Arc,
    },
    time::Duration,
};

#[derive(Clone)]
pub struct Timing {
    total_time_ms: Arc<AtomicUsize>,
    num_samples: Arc<AtomicUsize>,
}

impl Timing {
    pub fn new(duration: Duration) -> Self {
        let time_ms = Self::millis_u8(duration);
        Self {
            total_time_ms: Arc::new(AtomicUsize::new(time_ms)),
            num_samples: Arc::new(AtomicUsize::new(time_ms)),
        }
    }

    pub fn get_ms(&self) -> f64 {
        self.total_time_ms.load(Ordering::SeqCst) as f64
    }

    pub fn get_seconds(&self) -> f64 {
        self.total_time_ms.load(Ordering::SeqCst) as f64 / 1000.0
    }

    pub fn update(&mut self, duration: Duration) {
        let time_ms = Self::millis_u8(duration);
        self.total_time_ms.fetch_add(time_ms, Ordering::SeqCst);
        self.num_samples.fetch_add(time_ms, Ordering::SeqCst);
    }

    pub fn millis_u8(duration: Duration) -> usize {
        duration.as_millis().min(u64::MAX as u128) as usize
    }
}

impl Display for Timing {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{:.3}s", self.get_seconds())
    }
}

impl Default for Timing {
    fn default() -> Self {
        Self {
            total_time_ms: Arc::new(AtomicUsize::new(0)),
            num_samples: Arc::new(AtomicUsize::new(0)),
        }
    }
}

#[derive(Default, Clone)]
pub struct Timings {
    input: Timing,
    output: Timing,
    mmseqs: Timing,
    hmm_build: Timing,
    cloud_search: Timing,
    forward: Timing,
}

impl Timings {
    pub fn total_seconds(&self) -> f64 {
        self.input.get_seconds()
            + self.output.get_seconds()
            + self.mmseqs.get_seconds()
            + self.hmm_build.get_seconds()
            + self.cloud_search.get_seconds()
            + self.forward.get_seconds()
    }
}

#[derive(Clone, Default)]
pub struct Stats {
    num_samples: Arc<AtomicUsize>,
    num_passed_cloud: Arc<AtomicUsize>,
    num_passed_forward: Arc<AtomicUsize>,
    timings: Timings,
}

impl Stats {
    pub fn increment_samples(&mut self) {
        self.num_samples.fetch_add(1, Ordering::SeqCst);
    }

    pub fn increment_passed_cloud(&mut self) {
        self.num_samples.fetch_add(1, Ordering::SeqCst);
    }

    pub fn increment_passed_forward(&mut self) {
        self.num_samples.fetch_add(1, Ordering::SeqCst);
    }

    pub fn set_mmseqs_time(&mut self, time: Duration) {
        self.timings.mmseqs = Timing::new(time);
    }

    pub fn add_input_time(&mut self, time: Duration) {
        self.timings.input.update(time);
    }

    pub fn add_output_time(&mut self, time: Duration) {
        self.timings.output.update(time);
    }

    pub fn add_hmm_build_time(&mut self, time: Duration) {
        self.timings.hmm_build.update(time);
    }

    pub fn add_cloud_search_time(&mut self, time: Duration) {
        self.timings.cloud_search.update(time);
    }

    pub fn add_forward_time(&mut self, time: Duration) {
        self.timings.forward.update(time);
    }

    pub fn write(&mut self, out: &mut impl Write) -> anyhow::Result<()> {
        let total = self.timings.total_seconds();
        let input = self.timings.input.get_seconds();
        let output = self.timings.output.get_seconds();
        let mmseqs = self.timings.mmseqs.get_seconds();
        let cloud = self.timings.cloud_search.get_seconds();
        let forward = self.timings.forward.get_seconds();

        let width = format!("{total:.2}").len();
        writeln!(
            out,
            "input:  {:w$.2}s\t{:5.2}%",
            input,
            input / total * 100.0,
            w = width
        )?;
        writeln!(
            out,
            "output: {:w$.2}s\t{:5.2}%",
            output,
            output / total * 100.0,
            w = width
        )?;
        writeln!(
            out,
            "mmseqs: {:w$.2}s\t{:5.2}%",
            mmseqs,
            mmseqs / total * 100.0,
            w = width
        )?;
        writeln!(
            out,
            "cloud:  {:w$.2}s\t{:5.2}%",
            cloud,
            cloud / total * 100.0,
            w = width
        )?;
        writeln!(
            out,
            "fwd:    {:w$.2}s\t{:5.2}%",
            forward,
            forward / total * 100.0,
            w = width
        )?;
        writeln!(out, "total:  {:w$.2}s", total, w = width)?;

        Ok(())
    }
}
