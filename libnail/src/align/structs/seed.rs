use serde::{Deserialize, Serialize};

#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct Seed {
    pub seq_start: usize,
    pub seq_end: usize,
    pub prf_start: usize,
    pub prf_end: usize,
    pub score: f32,
}
