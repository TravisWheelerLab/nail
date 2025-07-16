use std::fmt::Display;

use serde::{Deserialize, Serialize};

#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct Seed {
    pub seq_start: usize,
    pub seq_end: usize,
    pub prf_start: usize,
    pub prf_end: usize,
    pub score: f32,
}

impl Display for Seed {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}:{}:{}:{}:{}",
            self.seq_start, self.seq_end, self.prf_start, self.prf_end, self.score
        )
    }
}
