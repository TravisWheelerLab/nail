pub mod constants {
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
    states: Vec<usize>,
    /// node index
    k: Vec<usize>,
    /// position emitted in seq
    i: Vec<usize>,
    posterior_probabilities: Vec<f32>,
}

impl Trace {
    pub fn new(profile_length: usize, sequence_length: usize) -> Self {
        Trace {
            length: 0,
            profile_length,
            sequence_length,
            states: vec![],
            k: vec![],
            i: vec![],
            posterior_probabilities: vec![],
        }
    }

    pub fn append_with_posterior_probability(
        &mut self,
        state: usize,
        k: usize,
        i: usize,
        posterior_probability: f32,
    ) {
        match state {
            // TRACE_N => {}
            // TRACE_C => {}
            TRACE_J => {
                if self.states[self.length - 1] == state {
                    self.i[self.length] = i;
                    self.posterior_probabilities[self.length] = posterior_probability;
                } else {
                    self.i[self.length] = 0;
                    self.posterior_probabilities[self.length] = 0.0;
                }
                self.k[self.length] = 0;
            }
            // TRACE_X => {}
            // TRACE_S => {}
            // TRACE_B => {}
            // TRACE_E => {}
            TRACE_T => {
                self.i[self.length] = 0;
                self.posterior_probabilities[self.length] = 0.0;
                self.k[self.length] = 0;
            }
            TRACE_D => {
                self.i[self.length] = 0;
                self.posterior_probabilities[self.length] = 0.0;
                self.k[self.length] = k;
            }
            // TRACE_M => {}
            TRACE_I => {
                self.i[self.length] = i;
                self.posterior_probabilities[self.length] = posterior_probability;
                self.k[self.length] = k;
            }
            _ => {
                panic!("no such state");
            }
        }
        self.states[self.length] = state;
        self.length += 1;
    }
}
