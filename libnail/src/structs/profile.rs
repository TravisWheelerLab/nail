use crate::alphabet::{
    AMINO_ALPHABET_WITH_DEGENERATE, AMINO_BACKGROUND_FREQUENCIES, AMINO_INVERSE_MAP,
    AMINO_INVERSE_MAP_LOWER, UTF8_SPACE,
};
use crate::structs::hmm::constants::{
    HMM_DELETE_TO_DELETE, HMM_DELETE_TO_MATCH, HMM_INSERT_TO_INSERT, HMM_INSERT_TO_MATCH,
    HMM_MATCH_TO_DELETE, HMM_MATCH_TO_INSERT, HMM_MATCH_TO_MATCH,
};
use crate::structs::hmm::Alphabet;
use crate::structs::trace::constants::{
    TRACE_B, TRACE_C, TRACE_D, TRACE_E, TRACE_I, TRACE_J, TRACE_M, TRACE_N, TRACE_S, TRACE_T,
};
use crate::structs::Hmm;
use crate::util::{f32_vec_argmax, LogAbuse};

use std::fmt;
use std::fmt::Formatter;

#[derive(Clone)]
pub struct Profile {
    /// The name of the profile
    pub name: String,
    /// The accession number of the profile
    pub accession: String,
    /// Model length (number of nodes)
    pub length: usize,
    /// Current target sequence length
    pub target_length: usize,
    /// Calculated upper bound on max sequence length
    pub max_length: usize,
    /// Transition scores
    pub transitions: Vec<[f32; 8]>,
    /// Match scores
    pub match_scores: Vec<Vec<f32>>,
    /// Insert scores
    pub insert_scores: Vec<Vec<f32>>,
    /// Transitions from special states (E, N, B, J, C)
    pub special_transitions: [[f32; 2]; 5],
    /// The expected number of times that the J state is used
    pub expected_j_uses: f32,
    /// The profile's consensus sequence
    pub consensus_sequence: Vec<u8>,
    /// The sequence alphabet
    pub alphabet: Alphabet,
    pub forward_tau: f32,
    pub forward_lambda: f32,
}

impl Profile {
    pub const LN_2: f32 = 0.69314718055994529;
    pub const LN_2_R: f32 = 1.44269504088896341;

    pub const MAX_ALPHABET_SIZE: usize = 20;
    pub const MAX_DEGENERATE_ALPHABET_SIZE: usize = 29;
    pub const GAP_INDEX: usize = 20;
    // non-residue character: "*"
    pub const NON_RESIDUE_IDX: usize = 27;
    // missing data character: "~"
    pub const MISSING_DATA_IDX: usize = 28;

    // special state indices
    pub const NUM_SPECIAL_STATES: usize = 5;
    pub const SPECIAL_E_IDX: usize = 0;
    pub const SPECIAL_N_IDX: usize = 1;
    pub const SPECIAL_J_IDX: usize = 2;
    pub const SPECIAL_B_IDX: usize = 3;
    pub const SPECIAL_C_IDX: usize = 4;

    pub const SPECIAL_STATE_IDX_TO_NAME: [&'static str; 5] = ["E", "N", "J", "B", "C"];

    // special transition indices
    pub const SPECIAL_LOOP_IDX: usize = 0;
    pub const SPECIAL_MOVE_IDX: usize = 1;

    /// The number of allowed state transitions under the model.
    pub const NUM_STATE_TRANSITIONS: usize = 8;
    pub const MATCH_TO_MATCH_IDX: usize = 0;
    pub const INSERT_TO_MATCH_IDX: usize = 1;
    pub const DELETE_TO_MATCH_IDX: usize = 2;
    pub const BEGIN_TO_MATCH_IDX: usize = 3;
    pub const MATCH_TO_DELETE_IDX: usize = 4;
    pub const DELETE_TO_DELETE_IDX: usize = 5;
    pub const MATCH_TO_INSERT_IDX: usize = 6;
    pub const INSERT_TO_INSERT_IDX: usize = 7;

    pub fn new(hmm: &Hmm) -> Self {
        let mut profile = Profile {
            name: hmm.header.name.clone(),
            accession: hmm.header.accession_number.clone(),
            length: hmm.header.model_length,
            target_length: 0,
            max_length: 0,
            transitions: vec![[0.0; 8]; hmm.header.model_length + 1],
            match_scores: vec![
                vec![0.0; Profile::MAX_DEGENERATE_ALPHABET_SIZE];
                hmm.header.model_length + 1
            ],
            insert_scores: vec![
                vec![0.0; Profile::MAX_DEGENERATE_ALPHABET_SIZE];
                hmm.header.model_length + 1
            ],
            special_transitions: [[0.0; 2]; 5],
            expected_j_uses: 0.0,
            // buffered with a space so that indexing starts at 1
            consensus_sequence: vec![UTF8_SPACE],
            alphabet: Alphabet::Amino,
            forward_tau: hmm.stats.forward_tau,
            forward_lambda: hmm.stats.forward_lambda,
        };

        for state in 0..Profile::NUM_STATE_TRANSITIONS {
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
            profile.transitions[model_position_idx - 1][Profile::BEGIN_TO_MATCH_IDX] =
                (match_occupancy[model_position_idx] / z).ln();
        }

        // these settings are for the non-multi-hit mode
        // N, C, and J transitions are set later by length config
        profile.special_transitions[Profile::SPECIAL_E_IDX][Profile::SPECIAL_MOVE_IDX] = 0.0;
        profile.special_transitions[Profile::SPECIAL_E_IDX][Profile::SPECIAL_LOOP_IDX] =
            -f32::INFINITY;
        profile.expected_j_uses = 0.0;

        // transition scores
        for i in 1..=profile.length {
            profile.transitions[i][Profile::MATCH_TO_MATCH_IDX] =
                hmm.model.transition_probabilities[i][HMM_MATCH_TO_MATCH].ln_or_inf();
            profile.transitions[i][Profile::MATCH_TO_INSERT_IDX] =
                hmm.model.transition_probabilities[i][HMM_MATCH_TO_INSERT].ln_or_inf();
            profile.transitions[i][Profile::MATCH_TO_DELETE_IDX] =
                hmm.model.transition_probabilities[i][HMM_MATCH_TO_DELETE].ln_or_inf();
            profile.transitions[i][Profile::INSERT_TO_MATCH_IDX] =
                hmm.model.transition_probabilities[i][HMM_INSERT_TO_MATCH].ln_or_inf();
            profile.transitions[i][Profile::INSERT_TO_INSERT_IDX] =
                hmm.model.transition_probabilities[i][HMM_INSERT_TO_INSERT].ln_or_inf();
            profile.transitions[i][Profile::DELETE_TO_MATCH_IDX] =
                hmm.model.transition_probabilities[i][HMM_DELETE_TO_MATCH].ln_or_inf();
            profile.transitions[i][Profile::DELETE_TO_DELETE_IDX] =
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

            for alphabet_idx in 0..Profile::MAX_ALPHABET_SIZE {
                // score is match ln(emission / background)
                // TODO: probably should make these casts unnecessary
                profile.match_scores[model_position_idx][alphabet_idx] =
                    (hmm.model.match_probabilities[model_position_idx][alphabet_idx] as f64
                        / AMINO_BACKGROUND_FREQUENCIES[alphabet_idx] as f64)
                        .ln() as f32;
            }
            // for the rest of the alphabet, we don't have scores from the HMM file
            profile.match_scores[model_position_idx][Profile::GAP_INDEX] = -f32::INFINITY;
            profile.match_scores[model_position_idx][Profile::NON_RESIDUE_IDX] = -f32::INFINITY;
            profile.match_scores[model_position_idx][Profile::MISSING_DATA_IDX] = -f32::INFINITY;

            // set the the rest of the degenerate characters
            for alphabet_idx in
                Profile::MAX_ALPHABET_SIZE..Profile::MAX_DEGENERATE_ALPHABET_SIZE - 3
            {
                let mut result: f32 = 0.0;
                let mut denominator: f32 = 0.0;
                for i in 0..Profile::MAX_ALPHABET_SIZE {
                    result += profile.match_scores[model_position_idx][i]
                        * AMINO_BACKGROUND_FREQUENCIES[i];
                    denominator += AMINO_BACKGROUND_FREQUENCIES[i];
                }
                profile.match_scores[model_position_idx][alphabet_idx] = result / denominator;
            }
        }

        // insert scores
        for alphabet_idx in 0..Profile::MAX_DEGENERATE_ALPHABET_SIZE {
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
            profile.insert_scores[model_position_idx][Profile::GAP_INDEX] = -f32::INFINITY;
            profile.insert_scores[model_position_idx][Profile::NON_RESIDUE_IDX] = -f32::INFINITY;
            profile.insert_scores[model_position_idx][Profile::MISSING_DATA_IDX] = -f32::INFINITY;
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
                TRACE_B => {
                    self.special_transition_score(Profile::SPECIAL_N_IDX, Profile::SPECIAL_MOVE_IDX)
                }
                TRACE_N => {
                    self.special_transition_score(Profile::SPECIAL_N_IDX, Profile::SPECIAL_LOOP_IDX)
                }
                _ => panic!(),
            },
            TRACE_B => match state_to {
                TRACE_M => self.transition_score(Profile::BEGIN_TO_MATCH_IDX, idx_to - 1),
                _ => panic!(),
            },
            TRACE_M => match state_to {
                TRACE_M => self.transition_score(Profile::MATCH_TO_MATCH_IDX, idx_from),
                TRACE_I => self.transition_score(Profile::MATCH_TO_INSERT_IDX, idx_from),
                TRACE_D => self.transition_score(Profile::MATCH_TO_DELETE_IDX, idx_from),
                TRACE_E => 0.0,
                _ => panic!(),
            },
            TRACE_D => match state_to {
                TRACE_M => self.transition_score(Profile::DELETE_TO_MATCH_IDX, idx_from),
                TRACE_D => self.transition_score(Profile::DELETE_TO_DELETE_IDX, idx_from),
                TRACE_E => 0.0,
                _ => panic!(),
            },
            TRACE_I => match state_to {
                TRACE_M => self.transition_score(Profile::INSERT_TO_MATCH_IDX, idx_from),
                TRACE_I => self.transition_score(Profile::INSERT_TO_INSERT_IDX, idx_from),
                _ => panic!(),
            },
            TRACE_E => match state_to {
                TRACE_C => {
                    self.special_transition_score(Profile::SPECIAL_E_IDX, Profile::SPECIAL_MOVE_IDX)
                }
                TRACE_J => {
                    self.special_transition_score(Profile::SPECIAL_E_IDX, Profile::SPECIAL_LOOP_IDX)
                }
                _ => panic!(),
            },
            TRACE_J => match state_to {
                TRACE_B => {
                    self.special_transition_score(Profile::SPECIAL_J_IDX, Profile::SPECIAL_MOVE_IDX)
                }
                TRACE_J => {
                    self.special_transition_score(Profile::SPECIAL_J_IDX, Profile::SPECIAL_LOOP_IDX)
                }
                _ => panic!(),
            },
            TRACE_C => match state_to {
                TRACE_T => {
                    self.special_transition_score(Profile::SPECIAL_C_IDX, Profile::SPECIAL_MOVE_IDX)
                }
                TRACE_C => {
                    self.special_transition_score(Profile::SPECIAL_C_IDX, Profile::SPECIAL_LOOP_IDX)
                }
                _ => panic!(),
            },
            _ => panic!(),
        }
    }

    /// Sets the length of the current target sequence to which the profile will be aligned.
    ///
    /// This also adjusts the loop and move transition scores for the special states N, J, C.
    pub fn configure_for_target_length(&mut self, length: usize) {
        self.target_length = length;

        let move_probability: f32 =
            (2.0 + self.expected_j_uses) / (length as f32 + 2.0 + self.expected_j_uses);

        let loop_probability: f32 = 1.0 - move_probability;

        let loop_score = loop_probability.ln();
        let move_score = move_probability.ln();

        self.special_transitions[Profile::SPECIAL_N_IDX][Profile::SPECIAL_LOOP_IDX] = loop_score;
        self.special_transitions[Profile::SPECIAL_J_IDX][Profile::SPECIAL_LOOP_IDX] = loop_score;
        self.special_transitions[Profile::SPECIAL_C_IDX][Profile::SPECIAL_LOOP_IDX] = loop_score;

        self.special_transitions[Profile::SPECIAL_N_IDX][Profile::SPECIAL_MOVE_IDX] = move_score;
        self.special_transitions[Profile::SPECIAL_J_IDX][Profile::SPECIAL_MOVE_IDX] = move_score;
        self.special_transitions[Profile::SPECIAL_C_IDX][Profile::SPECIAL_MOVE_IDX] = move_score;
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

            for _ in 0..Profile::MAX_DEGENERATE_ALPHABET_SIZE {
                write!(f, "  ----- ")?;
            }
            writeln!(f)?;

            for j in 0..Profile::MAX_DEGENERATE_ALPHABET_SIZE {
                write!(f, "{:8.4} ", self.match_scores[i][j])?;
            }
            writeln!(f)?;

            for j in 0..Profile::MAX_DEGENERATE_ALPHABET_SIZE {
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
