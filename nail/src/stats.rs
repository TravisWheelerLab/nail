use std::time::Duration;

#[derive(Clone)]
pub struct Timing {
    num_samples: u32,
    total: Duration,
    mean: Duration,
    min: Duration,
    max: Duration,
}

impl Timing {
    pub fn new(time: Duration) -> Self {
        Self {
            num_samples: 1,
            total: time,
            mean: time,
            min: time,
            max: time,
        }
    }

    pub fn update(&mut self, time: Duration) {
        self.total += time;
        self.min = self.min.min(time);
        self.max = self.max.max(time);
        self.mean = ((self.mean * self.num_samples) + time) / (self.num_samples + 1);
        self.num_samples += 1;
    }

    pub fn merge(&mut self, other: &Self) -> Self {
        let num_samples = self.num_samples + other.num_samples;
        Self {
            num_samples,
            total: self.total + other.total,
            mean: (self.total * self.num_samples + other.total * other.num_samples) / num_samples,
            min: self.min.min(other.min),
            max: self.max.max(other.max),
        }
    }
}

impl Default for Timing {
    fn default() -> Self {
        Self {
            num_samples: 0,
            total: Duration::new(0, 0),
            mean: Duration::new(0, 0),
            min: Duration::new(u64::MAX, u32::MAX),
            max: Duration::new(0, 0),
        }
    }
}

#[derive(Default)]
pub struct Timings {
    mmseqs: Timing,
    hmm_build: Timing,
    cloud_search: Timing,
    forward: Timing,
    traceback: Timing,
}

impl Timings {
    pub fn merge(&mut self, other: &Self) -> Self {
        Self {
            mmseqs: self.mmseqs.merge(&other.mmseqs),
            hmm_build: self.hmm_build.merge(&other.hmm_build),
            cloud_search: self.cloud_search.merge(&other.cloud_search),
            forward: self.forward.merge(&other.forward),
            traceback: self.traceback.merge(&other.traceback),
        }
    }
}

#[derive(Default)]
pub struct Stats {
    num_samples: usize,
    num_passed_cloud: usize,
    num_passed_forward: usize,
    timings: Timings,
}

impl Stats {
    pub fn merge(&mut self, other: &Stats) -> Self {
        Stats {
            num_samples: self.num_samples + other.num_samples,
            num_passed_cloud: self.num_passed_cloud + other.num_passed_cloud,
            num_passed_forward: self.num_passed_forward + other.num_passed_forward,
            timings: self.timings.merge(&other.timings),
        }
    }

    pub fn increment_passed_cloud(&mut self) {
        self.num_passed_cloud += 1;
    }

    pub fn increment_passed_forward(&mut self) {
        self.num_passed_forward += 1;
    }

    pub fn set_mmseqs_time(&mut self, time: Duration) {
        self.timings.mmseqs = Timing::new(time);
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

    pub fn add_traceback_time(&mut self, time: Duration) {
        self.timings.traceback.update(time);
    }
}
