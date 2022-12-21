mod constants {
    pub const TRACE_BOGUS: usize = 0;
    pub const TRACE_M: usize = 1;
    pub const TRACE_D: usize = 2;
    pub const TRACE_I: usize = 3;
    pub const TRACE_S: usize = 4;
    pub const TRACE_N: usize = 5;
    pub const TRACE_B: usize = 6;
    pub const TRACE_E: usize = 7;
    pub const TRACE_C: usize = 8;
    pub const TRACE_T: usize = 9;
    pub const TRACE_J: usize = 10;
    pub const TRACE_X: usize = 11;
}

use constants::*;

pub struct Trace {
    length: usize,
    profile_length: usize,
    sequence_length: usize,
    states: Vec<u8>,
    k: Vec<usize>,
    i: Vec<usize>,
    posterior_probabilities: Vec<f32>,
}

impl Trace {}
