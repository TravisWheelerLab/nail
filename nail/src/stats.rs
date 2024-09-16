use std::{
    io::Write,
    sync::{
        atomic::{AtomicU64, Ordering},
        Arc,
    },
    time::Duration,
};

#[repr(usize)]
#[derive(Clone, Copy)]
pub enum SerialTimed {
    Total,
    // Setup,
    Seeding,
    Alignment,
    _SENTINEL,
}
const NUM_SERIAL_TIMED: usize = SerialTimed::_SENTINEL as usize;

#[repr(usize)]
#[derive(Clone, Copy)]
pub enum ThreadTimed {
    Total,
    // Input,
    Output,
    MemoryInit,
    HmmBuild,
    CloudSearch,
    // BoundsAdjustment,
    Forward,
    Backward,
    Posterior,
    Traceback,
    NullTwo,
    _SENTINEL,
}
const NUM_THREAD_TIMED: usize = ThreadTimed::_SENTINEL as usize;

#[repr(usize)]
#[derive(Clone, Copy)]
pub enum FilterStage {
    Seed,
    Cloud,
    Forward,
    _SENTINEL,
}
const NUM_FILTERS: usize = FilterStage::_SENTINEL as usize;

#[derive(Clone, Default)]
pub struct Stats {
    serial_times: [Duration; NUM_SERIAL_TIMED],
    threaded_times: Arc<[AtomicU64; NUM_THREAD_TIMED]>,
    threaded_counts: Arc<[AtomicU64; NUM_THREAD_TIMED]>,
    filter_counts: Arc<[AtomicU64; NUM_FILTERS]>,
}

impl Stats {
    pub fn set_serial_time(&mut self, step: SerialTimed, time: Duration) {
        self.serial_times[step as usize] = time;
    }

    pub fn add_threaded_time(&mut self, step: ThreadTimed, time: Duration) {
        let time_nanos = Self::nanos(time);
        self.threaded_times[step as usize].fetch_add(time_nanos, Ordering::SeqCst);
        self.threaded_counts[step as usize].fetch_add(1, Ordering::SeqCst);
    }

    pub fn serial_time_total(&self, step: SerialTimed) -> Duration {
        self.serial_times[step as usize]
    }

    pub fn thread_time_total(&self, step: ThreadTimed) -> Duration {
        match step {
            ThreadTimed::_SENTINEL => {
                let total = self.thread_time_total(ThreadTimed::Total);

                let sum = self.threaded_times[1..]
                    .iter()
                    .map(|t| Duration::from_nanos(t.load(Ordering::SeqCst)))
                    .sum();

                total - sum
            }
            _ => {
                let nanos = self.threaded_times[step as usize].load(Ordering::SeqCst);
                let secs = nanos / 1e9 as u64;
                let nanos = (nanos - (secs * 1e9 as u64)) as u32;
                Duration::new(secs, nanos)
            }
        }
    }

    pub fn serial_time_pct(&self, step: SerialTimed) -> f64 {
        let total_nanos = Self::nanos(self.serial_times[ThreadTimed::Total as usize]) as f64;
        let nanos = Self::nanos(self.serial_times[step as usize]) as f64;

        nanos / total_nanos
    }

    pub fn thread_time_pct(&self, step: ThreadTimed) -> f64 {
        let total_nanos =
            self.threaded_times[ThreadTimed::Total as usize].load(Ordering::SeqCst) as f64;
        let nanos = match step {
            ThreadTimed::_SENTINEL => Self::nanos(self.thread_time_total(ThreadTimed::_SENTINEL)),
            _ => self.threaded_times[step as usize].load(Ordering::SeqCst),
        } as f64;

        nanos / total_nanos
    }

    pub fn increment_passed(&mut self, stage: FilterStage) {
        self.filter_counts[stage as usize].fetch_add(1, Ordering::SeqCst);
    }

    pub fn serial_string(&self, step: SerialTimed) -> String {
        let width = format!(
            "{:.2}",
            self.serial_time_total(SerialTimed::Total).as_secs_f64()
        )
        .len();

        format!(
            "{:w$.2}s ({:>5.2}%)",
            self.serial_time_total(step).as_secs_f64(),
            self.serial_time_pct(step) * 100.0,
            w = width,
        )
    }

    pub fn thread_string(&self, step: ThreadTimed) -> String {
        let width = format!(
            "{:.2}",
            self.thread_time_total(ThreadTimed::Total).as_secs_f64()
        )
        .len();

        format!(
            "{:w$.2}s ({:>5.2}%)",
            self.thread_time_total(step).as_secs_f64(),
            self.thread_time_pct(step) * 100.0,
            w = width,
        )
    }

    pub fn write(&self, out: &mut impl Write) -> anyhow::Result<()> {
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

        // writeln!(
        //     out,
        //     "     ├─ thread total:  {}",
        //     self.thread_string(ThreadTimed::Total)
        // )?;

        writeln!(
            out,
            "     ├─ memory init:  {}",
            self.thread_string(ThreadTimed::MemoryInit)
        )?;

        if self.thread_time_total(ThreadTimed::HmmBuild).as_secs_f64() > 0.0 {
            writeln!(
                out,
                "     ├─ hmm build:    {}",
                self.thread_string(ThreadTimed::HmmBuild)
            )?;
        }

        writeln!(
            out,
            "     ├─ cloud search: {}",
            self.thread_string(ThreadTimed::CloudSearch)
        )?;

        writeln!(
            out,
            "     ├─ forward:      {}",
            self.thread_string(ThreadTimed::Forward)
        )?;

        writeln!(
            out,
            "     ├─ backward:     {}",
            self.thread_string(ThreadTimed::Backward)
        )?;

        writeln!(
            out,
            "     ├─ posterior:    {}",
            self.thread_string(ThreadTimed::Posterior)
        )?;

        writeln!(
            out,
            "     ├─ traceback:    {}",
            self.thread_string(ThreadTimed::Traceback)
        )?;

        writeln!(
            out,
            "     ├─ null two:     {}",
            self.thread_string(ThreadTimed::NullTwo)
        )?;

        writeln!(
            out,
            "     └─ output:       {}",
            self.thread_string(ThreadTimed::Output)
        )?;

        writeln!(
            out,
            "     └─ [???]         {}",
            self.thread_string(ThreadTimed::_SENTINEL)
        )?;

        Ok(())
    }

    pub fn nanos(time: Duration) -> u64 {
        // u64::MAX nanoseconds is like 5,000,000 hours
        // or something, so this clamp should be fine.
        time.as_nanos().min(u64::MAX as u128) as u64
    }
}
