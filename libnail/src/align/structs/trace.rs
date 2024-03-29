use crate::alphabet::AMINO_INVERSE_MAP;
use crate::structs::{Profile, Sequence};
use std::io::Write;

#[derive(Default)]
pub struct Trace {
    pub length: usize,
    pub profile_length: usize,
    pub target_length: usize,
    pub states: Vec<usize>,
    /// node index
    pub profile_idx: Vec<usize>,
    /// position emitted in seq
    pub target_idx: Vec<usize>,
    pub posterior_probabilities: Vec<f32>,
    ///
    pub cell_fraction: f32,
}

impl Trace {
    pub const TRACE_IDX_TO_NAME: [&str; 12] = [
        "INVALID", "M", "D", "I", "S", "N", "B", "E", "C", "T", "J", "X",
    ];

    pub const INVALID_STATE: usize = 0;
    pub const M_STATE: usize = 1;
    pub const D_STATE: usize = 2;
    pub const I_STATE: usize = 3;
    pub const S_STATE: usize = 4;
    pub const N_STATE: usize = 5;
    pub const B_STATE: usize = 6;
    pub const E_STATE: usize = 7;
    pub const C_STATE: usize = 8;
    pub const T_STATE: usize = 9;
    pub const J_STATE: usize = 10;
    pub const X_STATE: usize = 11;

    pub fn new(target_length: usize, profile_length: usize) -> Self {
        Trace {
            length: 0,
            profile_length,
            target_length,
            states: vec![],
            profile_idx: vec![],
            target_idx: vec![],
            posterior_probabilities: vec![],
            cell_fraction: 0.0,
        }
    }

    pub fn resize(&mut self, new_target_length: usize, new_profile_length: usize) {
        #![allow(unused_variables)]
        todo!();
    }

    pub fn append_with_posterior_probability(
        &mut self,
        state: usize,
        target_idx: usize,
        profile_idx: usize,
        posterior_probability: f32,
    ) {
        match state {
            Trace::N_STATE | Trace::C_STATE | Trace::J_STATE => {
                if self.states[self.length - 1] == state {
                    self.target_idx.push(target_idx);
                    self.posterior_probabilities.push(posterior_probability);
                } else {
                    self.target_idx.push(0);
                    self.posterior_probabilities.push(0.0);
                }
                self.profile_idx.push(0);
            }
            Trace::X_STATE | Trace::S_STATE | Trace::B_STATE | Trace::E_STATE | Trace::T_STATE => {
                self.target_idx.push(0);
                self.posterior_probabilities.push(0.0);
                self.profile_idx.push(0);
            }
            Trace::D_STATE => {
                self.target_idx.push(0);
                self.posterior_probabilities.push(0.0);
                self.profile_idx.push(profile_idx);
            }
            Trace::M_STATE | Trace::I_STATE => {
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
            if (self.states[z] == Trace::N_STATE && self.states[z + 1] == Trace::N_STATE)
                || (self.states[z] == Trace::C_STATE && self.states[z + 1] == Trace::C_STATE)
                || (self.states[z] == Trace::J_STATE && self.states[z + 1] == Trace::J_STATE)
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
                Trace::TRACE_IDX_TO_NAME[current_state],
                self.profile_idx[trace_idx],
                self.target_idx[trace_idx],
                transition_score
            )?;

            if current_state == Trace::M_STATE {
                write!(
                    out,
                    " {:8.4}",
                    profile.match_score(current_residue as usize, self.profile_idx[trace_idx])
                )?;

                bit_score +=
                    profile.match_score(current_residue as usize, self.profile_idx[trace_idx]);
                write!(out, " {:8.4}", self.posterior_probabilities[trace_idx])?;
                accuracy += self.posterior_probabilities[trace_idx];
            } else if current_state == Trace::I_STATE {
                write!(
                    out,
                    " {:8.4}",
                    profile.insert_score(current_residue as usize, self.profile_idx[trace_idx])
                )?;

                bit_score +=
                    profile.insert_score(current_residue as usize, self.profile_idx[trace_idx]);

                write!(out, " {:8.4}", self.posterior_probabilities[trace_idx])?;

                accuracy += self.posterior_probabilities[trace_idx];
            } else if (current_state == Trace::N_STATE
                && self.states[trace_idx - 1] == Trace::N_STATE)
                || (current_state == Trace::C_STATE && self.states[trace_idx - 1] == Trace::C_STATE)
                || (current_state == Trace::J_STATE && self.states[trace_idx - 1] == Trace::J_STATE)
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
