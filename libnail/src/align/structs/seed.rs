use serde::{Deserialize, Serialize};

#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct Seed {
    pub target_start: usize,
    pub target_end: usize,
    pub profile_start: usize,
    pub profile_end: usize,
}
