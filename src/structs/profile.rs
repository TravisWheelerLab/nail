use crate::alphabet::{
    AMINO_ALPHABET_WITH_DEGENERATE, AMINO_BACKGROUND_FREQUENCIES, AMINO_INVERSE_MAP,
    AMINO_INVERSE_MAP_LOWER, UTF8_SPACE,
};
use crate::structs::hmm::constants::{
    HMM_DELETE_TO_DELETE, HMM_DELETE_TO_MATCH, HMM_INSERT_TO_INSERT, HMM_INSERT_TO_MATCH,
    HMM_MATCH_TO_DELETE, HMM_MATCH_TO_INSERT, HMM_MATCH_TO_MATCH,
};
use crate::structs::hmm::P7Alphabet;
use crate::structs::trace::constants::{
    TRACE_B, TRACE_C, TRACE_D, TRACE_E, TRACE_I, TRACE_J, TRACE_M, TRACE_N, TRACE_S, TRACE_T,
};
use crate::structs::Hmm;
use crate::util::{f32_vec_argmax, LogAbuse, PrintMe};

use std::fmt;
use std::fmt::Formatter;
use std::io::Write;

use anyhow::Result;

pub mod constants {
    // this is p7_MAXABET; (4 or 20)
    pub const MAX_ALPHABET_SIZE: usize = 20;
    // this is p7_MAXCODE in hmmer; (18 or 29)
    pub const MAX_DEGENERATE_ALPHABET_SIZE: usize = 29;

    pub const GAP_INDEX: usize = 20;
    // non-residue character: appears to be "*"
    pub const NON_RESIDUE_INDEX: usize = 27;
    // missing data character: appears to be "~"
    pub const MISSING_DATA_INDEX: usize = 28;

    // TODO: maybe just use the std const?
    pub const LN_2: f32 = 0.69314718055994529;
    // const LN_2: f32 = std::f32::consts::LN_2;
    // const LN_2: f64 = std::f64::consts::LN_2;

    pub const LN_2_R: f32 = 1.44269504088896341;
    // const LN_2_R: f32 = std::f32::consts::LOG2_E;
    // const LN_2_R: f64 = std::f64::consts::LOG2_E;

    // special state indices
    pub const NUM_SPECIAL_STATES: usize = 5;
    pub const SPECIAL_E: usize = 0;
    pub const SPECIAL_N: usize = 1;
    pub const SPECIAL_J: usize = 2;
    pub const SPECIAL_B: usize = 3;
    pub const SPECIAL_C: usize = 4;

    pub const SPECIAL_STATE_IDX_TO_NAME: [&str; 5] = ["E", "N", "J", "B", "C"];

    // special transition indices
    pub const SPECIAL_LOOP: usize = 0;
    pub const SPECIAL_MOVE: usize = 1;

    // transition indices
    pub const PROFILE_NUM_TRANSITIONS: usize = 8;
    pub const PROFILE_MATCH_TO_MATCH: usize = 0;
    pub const PROFILE_INSERT_TO_MATCH: usize = 1;
    pub const PROFILE_DELETE_TO_MATCH: usize = 2;
    pub const PROFILE_BEGIN_TO_MATCH: usize = 3;
    pub const PROFILE_MATCH_TO_DELETE: usize = 4;
    pub const PROFILE_DELETE_TO_DELETE: usize = 5;
    pub const PROFILE_MATCH_TO_INSERT: usize = 6;
    pub const PROFILE_INSERT_TO_INSERT: usize = 7;
}
use constants::*;

pub struct Profile {
    /// The name of the profile
    pub name: String,
    /// M: model length (number of nodes)
    pub length: usize,
    /// T: current target sequence length
    pub target_length: usize,
    /// Calculated upper bound on max sequence length
    pub max_length: usize,
    /// tsc: Transition scores
    pub transitions: Vec<[f32; 8]>,
    /// rsc: Match scores
    pub match_scores: Vec<Vec<f32>>,
    /// rsc: Insert scores
    pub insert_scores: Vec<Vec<f32>>,
    /// xsc: Transitions from special states (E, N, B, J, C)
    pub special_transitions: [[f32; 2]; 5],
    /// The expected number of times that the J state is used
    pub expected_j_uses: f32,
    /// The profile's consensus sequence
    pub consensus_sequence: Vec<u8>,
    /// The sequence alphabet
    pub alphabet: P7Alphabet,
}

impl Profile {
    pub fn new(hmm: &Hmm) -> Self {
        let mut profile = Profile {
            name: hmm.header.name.clone(),
            length: hmm.header.model_length,
            target_length: 0,
            max_length: 0,
            transitions: vec![[0.0; 8]; hmm.header.model_length + 1],
            match_scores: vec![
                vec![0.0; MAX_DEGENERATE_ALPHABET_SIZE];
                hmm.header.model_length + 1
            ],
            insert_scores: vec![
                vec![0.0; MAX_DEGENERATE_ALPHABET_SIZE];
                hmm.header.model_length + 1
            ],
            special_transitions: [[0.0; 2]; 5],
            expected_j_uses: 0.0,
            consensus_sequence: vec![UTF8_SPACE],
            alphabet: P7Alphabet::Amino,
        };

        for state in 0..PROFILE_NUM_TRANSITIONS {
            profile.transitions[0][state] = -f32::INFINITY;
        }

        // porting p7_hmm_CalculateOccupancy() from p7_hmm.c
        let mut match_occupancy = vec![0.0; profile.length + 1];

        match_occupancy[1] = hmm.model.transition_probabilities[0][HMM_MATCH_TO_INSERT]
            + hmm.model.transition_probabilities[0][HMM_MATCH_TO_MATCH];

        for k in 2..=profile.length {
            match_occupancy[k] = match_occupancy[k - 1]
                * (hmm.model.transition_probabilities[k - 1][HMM_MATCH_TO_MATCH]
                    + hmm.model.transition_probabilities[k - 1][HMM_MATCH_TO_INSERT])
                + (1.0 - match_occupancy[k - 1])
                    * hmm.model.transition_probabilities[k - 1][HMM_DELETE_TO_MATCH]
        }

        // TODO: what does Z represent?
        let mut z: f32 = 0.0;

        for profile_idx in 1..=profile.length {
            z += match_occupancy[profile_idx] * (profile.length - profile_idx + 1) as f32;
        }

        // the goal here must be to set the transition of begin to match at each position
        for model_position_idx in 1..=profile.length {
            profile.transitions[model_position_idx - 1][PROFILE_BEGIN_TO_MATCH] =
                (match_occupancy[model_position_idx] / z).ln();
        }

        // these settings are for the non-multi-hit mode
        // N, C, and J transitions are set later by length config
        profile.special_transitions[SPECIAL_E][SPECIAL_MOVE] = 0.0;
        profile.special_transitions[SPECIAL_E][SPECIAL_LOOP] = -f32::INFINITY;
        profile.expected_j_uses = 0.0;

        // transition scores
        for i in 1..=profile.length {
            profile.transitions[i][PROFILE_MATCH_TO_MATCH] =
                hmm.model.transition_probabilities[i][HMM_MATCH_TO_MATCH].ln_or_inf();
            profile.transitions[i][PROFILE_MATCH_TO_INSERT] =
                hmm.model.transition_probabilities[i][HMM_MATCH_TO_INSERT].ln_or_inf();
            profile.transitions[i][PROFILE_MATCH_TO_DELETE] =
                hmm.model.transition_probabilities[i][HMM_MATCH_TO_DELETE].ln_or_inf();
            profile.transitions[i][PROFILE_INSERT_TO_MATCH] =
                hmm.model.transition_probabilities[i][HMM_INSERT_TO_MATCH].ln_or_inf();
            profile.transitions[i][PROFILE_INSERT_TO_INSERT] =
                hmm.model.transition_probabilities[i][HMM_INSERT_TO_INSERT].ln_or_inf();
            profile.transitions[i][PROFILE_DELETE_TO_MATCH] =
                hmm.model.transition_probabilities[i][HMM_DELETE_TO_MATCH].ln_or_inf();
            profile.transitions[i][PROFILE_DELETE_TO_DELETE] =
                hmm.model.transition_probabilities[i][HMM_DELETE_TO_DELETE].ln_or_inf();
        }

        // match scores
        for model_position_idx in 1..=profile.length {
            // the consensus residue is the match emission with the highest probability
            // TODO: this should not be the case for single sequence model
            // the argmax of this positions match probability vector will be the digital
            // index of the residue that has the highest match emission probability
            let match_probabilities_argmax: usize =
                f32_vec_argmax(&hmm.model.match_probabilities[model_position_idx]);
            let match_probabilities_max: f32 =
                hmm.model.match_probabilities[model_position_idx][match_probabilities_argmax];

            profile
                .consensus_sequence
                .push(if match_probabilities_max > 0.5 {
                    // if the match emission probability for the residue is greater
                    // than 0.50 (amino), we want to display it as a capital letter
                    AMINO_INVERSE_MAP[&(match_probabilities_argmax as u8)]
                } else {
                    // otherwise, we want to display it as a lowercase letter
                    AMINO_INVERSE_MAP_LOWER[&(match_probabilities_argmax as u8)]
                });

            for alphabet_idx in 0..MAX_ALPHABET_SIZE {
                // score is match ln(emission / background)
                // TODO: probably should make these casts unnecessary
                profile.match_scores[model_position_idx][alphabet_idx] =
                    (hmm.model.match_probabilities[model_position_idx][alphabet_idx] as f64
                        / AMINO_BACKGROUND_FREQUENCIES[alphabet_idx] as f64)
                        .ln() as f32;
            }
            // for the rest of the alphabet, we don't have scores from the HMM file
            profile.match_scores[model_position_idx][GAP_INDEX] = -f32::INFINITY;
            profile.match_scores[model_position_idx][NON_RESIDUE_INDEX] = -f32::INFINITY;
            profile.match_scores[model_position_idx][MISSING_DATA_INDEX] = -f32::INFINITY;

            // set the the rest of the degenerate characters
            for alphabet_idx in MAX_ALPHABET_SIZE..MAX_DEGENERATE_ALPHABET_SIZE - 3 {
                let mut result: f32 = 0.0;
                let mut denominator: f32 = 0.0;
                for i in 0..MAX_ALPHABET_SIZE {
                    result += profile.match_scores[model_position_idx][i]
                        * AMINO_BACKGROUND_FREQUENCIES[i];
                    denominator += AMINO_BACKGROUND_FREQUENCIES[i];
                }
                profile.match_scores[model_position_idx][alphabet_idx] = result / denominator;
            }
        }

        // insert scores
        for alphabet_idx in 0..MAX_DEGENERATE_ALPHABET_SIZE {
            for model_position_idx in 1..profile.length {
                // setting insert scores to 0 corresponds to insertion
                // emissions being equal to background probabilities
                //    ** because ln(P/P) = ln(1) = 0
                profile.insert_scores[model_position_idx][alphabet_idx] = 0.0;
            }
            // insert at position M should be impossible,
            profile.insert_scores[profile.length][alphabet_idx] = -f32::INFINITY;
        }

        for model_position_idx in 0..profile.length {
            profile.insert_scores[model_position_idx][GAP_INDEX] = -f32::INFINITY;
            profile.insert_scores[model_position_idx][NON_RESIDUE_INDEX] = -f32::INFINITY;
            profile.insert_scores[model_position_idx][MISSING_DATA_INDEX] = -f32::INFINITY;
        }

        profile
    }

    #[inline(always)]
    pub fn match_score(&self, alphabet_idx: usize, profile_idx: usize) -> f32 {
        self.match_scores[profile_idx][alphabet_idx]
    }

    #[inline(always)]
    pub fn insert_score(&self, alphabet_idx: usize, profile_idx: usize) -> f32 {
        self.insert_scores[profile_idx][alphabet_idx]
    }

    #[inline(always)]
    pub fn transition_score(&self, transition_idx: usize, profile_idx: usize) -> f32 {
        self.transitions[profile_idx][transition_idx]
    }

    /// This is essentially a Kronecker delta function that returns 1.0 when the transition
    /// score is finite and f32::MIN when the transition score is -f32::INFINITY.
    ///
    /// Semantically, this means we are disallowing "impossible state paths" during posterior traceback.
    #[inline(always)]
    pub fn transition_score_delta(&self, transition_idx: usize, profile_idx: usize) -> f32 {
        if self.transitions[profile_idx][transition_idx] == -f32::INFINITY {
            f32::MIN_POSITIVE
        } else {
            1.0
        }
    }

    #[inline(always)]
    pub fn special_transition_score(&self, state_idx: usize, transition_idx: usize) -> f32 {
        self.special_transitions[state_idx][transition_idx]
    }

    /// This is essentially a Kronecker delta function that returns 1.0 when the transition
    /// score is finite and f32::MIN when the transition score is -f32::INFINITY.
    ///
    /// Semantically, this means we are disallowing "impossible state paths" during posterior traceback.
    #[inline(always)]
    pub fn special_transition_score_delta(&self, state_idx: usize, transition_idx: usize) -> f32 {
        if self.special_transitions[state_idx][transition_idx] == -f32::INFINITY {
            f32::MIN_POSITIVE
        } else {
            1.0
        }
    }

    pub fn generic_transition_score(
        &self,
        state_from: usize,
        idx_from: usize,
        state_to: usize,
        idx_to: usize,
    ) -> f32 {
        match state_from {
            TRACE_S | TRACE_T => 0.0,
            TRACE_N => match state_to {
                TRACE_B => self.special_transition_score(SPECIAL_N, SPECIAL_MOVE),
                TRACE_N => self.special_transition_score(SPECIAL_N, SPECIAL_LOOP),
                _ => panic!(),
            },
            TRACE_B => match state_to {
                TRACE_M => self.transition_score(PROFILE_BEGIN_TO_MATCH, idx_to - 1),
                _ => panic!(),
            },
            TRACE_M => match state_to {
                TRACE_M => self.transition_score(PROFILE_MATCH_TO_MATCH, idx_from),
                TRACE_I => self.transition_score(PROFILE_MATCH_TO_INSERT, idx_from),
                TRACE_D => self.transition_score(PROFILE_MATCH_TO_DELETE, idx_from),
                TRACE_E => 0.0,
                _ => panic!(),
            },
            TRACE_D => match state_to {
                TRACE_M => self.transition_score(PROFILE_DELETE_TO_MATCH, idx_from),
                TRACE_D => self.transition_score(PROFILE_DELETE_TO_DELETE, idx_from),
                TRACE_E => 0.0,
                _ => panic!(),
            },
            TRACE_I => match state_to {
                TRACE_M => self.transition_score(PROFILE_INSERT_TO_MATCH, idx_from),
                TRACE_I => self.transition_score(PROFILE_INSERT_TO_INSERT, idx_from),
                _ => panic!(),
            },
            TRACE_E => match state_to {
                TRACE_C => self.special_transition_score(SPECIAL_E, SPECIAL_MOVE),
                TRACE_J => self.special_transition_score(SPECIAL_E, SPECIAL_LOOP),
                _ => panic!(),
            },
            TRACE_J => match state_to {
                TRACE_B => self.special_transition_score(SPECIAL_J, SPECIAL_MOVE),
                TRACE_J => self.special_transition_score(SPECIAL_J, SPECIAL_LOOP),
                _ => panic!(),
            },
            TRACE_C => match state_to {
                TRACE_T => self.special_transition_score(SPECIAL_C, SPECIAL_MOVE),
                TRACE_C => self.special_transition_score(SPECIAL_C, SPECIAL_LOOP),
                _ => panic!(),
            },
            _ => panic!(),
        }
    }

    /// Sets the length of the current target sequence to which the profile will be aligned.
    ///
    /// This also adjusts the loop and move transition scores for the special states N, J, C.
    pub fn configure_for_length(&mut self, length: usize) {
        // TODO: should we set up some sort of flag to try to make
        //       sure the model has been initialized for scoring?
        self.target_length = length;

        // a somewhat indecipherable comment from hmmer:
        //   2/(L+2) for sw; 3/(L+3) for fs
        let move_probability: f32 =
            (2.0 + self.expected_j_uses) / (length as f32 + 2.0 + self.expected_j_uses);

        let loop_probability: f32 = 1.0 - move_probability;

        let loop_score = loop_probability.ln();
        let move_score = move_probability.ln();

        self.special_transitions[SPECIAL_N][SPECIAL_LOOP] = loop_score;
        self.special_transitions[SPECIAL_J][SPECIAL_LOOP] = loop_score;
        self.special_transitions[SPECIAL_C][SPECIAL_LOOP] = loop_score;

        self.special_transitions[SPECIAL_N][SPECIAL_MOVE] = move_score;
        self.special_transitions[SPECIAL_J][SPECIAL_MOVE] = move_score;
        self.special_transitions[SPECIAL_C][SPECIAL_MOVE] = move_score;
    }
}

impl fmt::Debug for Profile {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        writeln!(f, "model length: {}", self.length)?;
        writeln!(f, "target length: {}", self.target_length)?;

        for i in 0..5 {
            writeln!(
                f,
                "{:8.4} {:8.4}",
                self.special_transitions[i][0], self.special_transitions[i][1]
            )?;
        }

        for i in 0..=self.length {
            writeln!(f, "{}", i)?;
            for residue in AMINO_ALPHABET_WITH_DEGENERATE {
                write!(f, "    {}   ", residue)?;
            }
            writeln!(f)?;

            for _ in 0..MAX_DEGENERATE_ALPHABET_SIZE {
                write!(f, "  ----- ")?;
            }
            writeln!(f)?;

            for j in 0..MAX_DEGENERATE_ALPHABET_SIZE {
                write!(f, "{:8.4} ", self.match_scores[i][j])?;
            }
            writeln!(f)?;

            for j in 0..MAX_DEGENERATE_ALPHABET_SIZE {
                write!(f, "{:8.4} ", self.insert_scores[i][j])?;
            }
            writeln!(f)?;

            for t in 0..8 {
                write!(f, "{:8.4} ", self.transitions[i][t])?;
            }
            writeln!(f)?;
            writeln!(f)?;
        }

        Ok(())
    }
}
