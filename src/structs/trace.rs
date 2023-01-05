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

    pub const TRACE_IDX_TO_NAME: [&str; 12] = [
        "BOGUS", "M", "D", "I", "S", "N", "B", "E", "C", "T", "J", "X",
    ];
}

use crate::alphabet::AMINO_INVERSE_MAP;
use crate::structs::{Profile, Sequence};
use constants::*;
use std::io::Write;

pub struct Trace {
    pub length: usize,
    pub profile_length: usize,
    pub sequence_length: usize,
    pub states: Vec<usize>,
    /// node index
    pub profile_idx: Vec<usize>,
    /// position emitted in seq
    pub target_idx: Vec<usize>,
    pub posterior_probabilities: Vec<f32>,
}

impl Trace {
    pub fn new(profile_length: usize, sequence_length: usize) -> Self {
        Trace {
            length: 0,
            profile_length,
            sequence_length,
            states: vec![],
            profile_idx: vec![],
            target_idx: vec![],
            posterior_probabilities: vec![],
        }
    }

    pub fn append_with_posterior_probability(
        &mut self,
        state: usize,
        profile_idx: usize,
        target_idx: usize,
        posterior_probability: f32,
    ) {
        match state {
            TRACE_N | TRACE_C | TRACE_J => {
                if self.states[self.length - 1] == state {
                    self.target_idx.push(target_idx);
                    self.posterior_probabilities.push(posterior_probability);
                } else {
                    self.target_idx.push(0);
                    self.posterior_probabilities.push(0.0);
                }
                self.profile_idx.push(0);
            }
            TRACE_X | TRACE_S | TRACE_B | TRACE_E | TRACE_T => {
                self.target_idx.push(0);
                self.posterior_probabilities.push(0.0);
                self.profile_idx.push(0);
            }
            TRACE_D => {
                self.target_idx.push(0);
                self.posterior_probabilities.push(0.0);
                self.profile_idx.push(profile_idx);
            }
            TRACE_M | TRACE_I => {
                self.target_idx.push(target_idx);
                self.posterior_probabilities.push(posterior_probability);
                self.profile_idx.push(profile_idx);
            }
            _ => {
                panic!("no such state {}", state);
            }
        }
        self.states.push(state);
        self.length += 1;
    }

    pub fn reverse(&mut self) {
        for z in 0..self.length {
            if (self.states[z] == TRACE_N && self.states[z + 1] == TRACE_N)
                || (self.states[z] == TRACE_C && self.states[z + 1] == TRACE_C)
                || (self.states[z] == TRACE_J && self.states[z + 1] == TRACE_J)
            {
                if self.target_idx[z] == 0 && self.target_idx[z + 1] > 0 {
                    self.target_idx[z] = self.target_idx[z + 1];
                    self.target_idx[z + 1] = 0;
                    self.posterior_probabilities[z] = self.posterior_probabilities[z + 1];
                    self.posterior_probabilities[z + 1] = 0.0;
                }
            }
        }

        self.states.reverse();
        self.target_idx.reverse();
        self.profile_idx.reverse();
        self.posterior_probabilities.reverse();
    }

    pub fn dump(
        &self,
        out: &mut impl Write,
        profile: &Profile,
        target: &Sequence,
    ) -> anyhow::Result<()> {
        let mut bit_score: f32 = 0.0;
        let mut accuracy: f32 = 0.0;

        writeln!(
            out,
            "st   p     t      transit emission postprob - traceback len {}",
            self.length
        )?;
        writeln!(out, "--  ---- ------  -------- -------- --------")?;
        for trace_idx in 0..self.length {
            let current_state = self.states[trace_idx];
            let current_residue = target.digital_bytes[self.target_idx[trace_idx]];

            let transition_score = if trace_idx < self.length - 1 {
                profile.generic_transition_score(
                    current_state,
                    self.profile_idx[trace_idx],
                    self.states[trace_idx + 1],
                    self.profile_idx[trace_idx + 1],
                )
            } else {
                0.0
            };
            bit_score += transition_score;

            write!(
                out,
                "{:1}  {:4} {:6}   {:8.4}",
                TRACE_IDX_TO_NAME[current_state],
                self.profile_idx[trace_idx],
                self.target_idx[trace_idx],
                transition_score
            )?;

            if current_state == TRACE_M {
                write!(
                    out,
                    " {:8.4}",
                    profile.match_score(current_residue as usize, self.profile_idx[trace_idx])
                )?;

                bit_score +=
                    profile.match_score(current_residue as usize, self.profile_idx[trace_idx]);
                write!(out, " {:8.4}", self.posterior_probabilities[trace_idx])?;
                accuracy += self.posterior_probabilities[trace_idx];
            } else if current_state == TRACE_I {
                write!(
                    out,
                    " {:8.4}",
                    profile.insert_score(current_residue as usize, self.profile_idx[trace_idx])
                )?;

                bit_score +=
                    profile.insert_score(current_residue as usize, self.profile_idx[trace_idx]);

                write!(out, " {:8.4}", self.posterior_probabilities[trace_idx])?;

                accuracy += self.posterior_probabilities[trace_idx];
            } else if (current_state == TRACE_N && self.states[trace_idx - 1] == TRACE_N)
                || (current_state == TRACE_C && self.states[trace_idx - 1] == TRACE_C)
                || (current_state == TRACE_J && self.states[trace_idx - 1] == TRACE_J)
            {
                write!(out, " {:8}", 0)?;
                write!(out, " {:8.4}", self.posterior_probabilities[trace_idx])?;
                accuracy += self.posterior_probabilities[trace_idx];
            }

            writeln!(out, " {}", AMINO_INVERSE_MAP[&current_residue] as char)?;
        }

        writeln!(out, "                -------- -------- --------")?;
        writeln!(
            out,
            "                  total: {:8.4} {:8.4}\n",
            bit_score, accuracy
        )?;

        Ok(())
    }
}
