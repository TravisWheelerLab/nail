use rand::SeedableRng;
use rand_pcg::Pcg64;

use crate::align::structs::{DpMatrixSparse, RowBounds, Trace};
use crate::align::{forward, null_one_score};
use crate::alphabet::{
    AMINO_ALPHABET_WITH_DEGENERATE, AMINO_BACKGROUND_FREQUENCIES, AMINO_INVERSE_MAP,
    AMINO_INVERSE_MAP_LOWER, UTF8_SPACE,
};
use crate::structs::hmm::constants::{
    HMM_DELETE_TO_DELETE, HMM_DELETE_TO_MATCH, HMM_INSERT_TO_INSERT, HMM_INSERT_TO_MATCH,
    HMM_MATCH_TO_DELETE, HMM_MATCH_TO_INSERT, HMM_MATCH_TO_MATCH,
};
use crate::structs::hmm::Alphabet;
use crate::structs::Hmm;
use crate::util::{avg_relative_entropy, f32_vec_argmax, LogAbuse, VecUtils};

use std::fmt;
use std::fmt::Formatter;

use super::Sequence;

#[derive(Clone, Default)]
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
    pub consensus_sequence_bytes_utf8: Vec<u8>,
    /// The sequence alphabet
    pub alphabet: Alphabet,
    pub forward_tau: f32,
    pub forward_lambda: f32,
}

impl Profile {
    // pub const LN_2: f32 = 0.69314718055994529;
    pub const LN_2: f32 = std::f32::consts::LN_2;
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

    pub fn tune_relative_entropy(&mut self, target: f32) {
        const TOLERANCE: f32 = 1e-3;

        let start_probs_by_pos: Vec<Vec<f32>> = self
            .match_scores
            .iter()
            .map(|scores| {
                scores
                    .iter()
                    .take(20)
                    .enumerate()
                    .map(|(idx, score)| score.exp() * AMINO_BACKGROUND_FREQUENCIES[idx])
                    .collect::<Vec<f32>>()
            })
            .collect();

        enum Mode {
            Growing,
            Binary,
        }

        let mut search_mode = Mode::Growing;

        let mut new_probs_by_pos = start_probs_by_pos.clone();
        let mut factor = 0.1;
        let mut upper_bound = 0.0;
        let mut lower_bound = 0.0;
        let mut relative_entropy =
            avg_relative_entropy(&new_probs_by_pos, &AMINO_BACKGROUND_FREQUENCIES);

        while (relative_entropy - target).abs() > TOLERANCE {
            new_probs_by_pos
                .iter_mut()
                .zip(&start_probs_by_pos)
                // skip model position 0
                .skip(1)
                .for_each(|(new_probs, start_probs)| {
                    new_probs
                        .iter_mut()
                        .zip(start_probs)
                        .enumerate()
                        .take(Profile::MAX_ALPHABET_SIZE)
                        .for_each(|(idx, (p_new, p_start))| {
                            *p_new = p_start - factor * AMINO_BACKGROUND_FREQUENCIES[idx]
                        });

                    new_probs.saturate_lower(1e-3);
                    new_probs.normalize();
                });

            relative_entropy =
                avg_relative_entropy(&new_probs_by_pos[1..], &AMINO_BACKGROUND_FREQUENCIES);

            match search_mode {
                Mode::Growing => {
                    if relative_entropy > target {
                        upper_bound = factor;
                        factor = (lower_bound + upper_bound) / 2.0;
                        search_mode = Mode::Binary;
                    } else {
                        lower_bound = factor;
                        factor += 0.1;
                    }
                }
                Mode::Binary => {
                    if relative_entropy > target {
                        upper_bound = factor;
                    } else {
                        lower_bound = factor;
                    }
                    factor = (lower_bound + upper_bound) / 2.0;
                }
            }
        }

        self.match_scores
            .iter_mut()
            .zip(new_probs_by_pos)
            // skip model position 0
            .skip(1)
            .for_each(|(scores, probs)| {
                scores
                    .iter_mut()
                    .zip(probs)
                    .enumerate()
                    .take(Profile::MAX_ALPHABET_SIZE)
                    .for_each(|(idx, (s, p))| *s = (p / AMINO_BACKGROUND_FREQUENCIES[idx]).ln())
            });
    }

    pub fn calibrate_tau(&mut self, n: usize, target_length: usize, tail_probability: f32) {
        self.configure_for_target_length(target_length);

        let mut row_bounds = RowBounds::new(target_length);
        row_bounds.fill_rectangle(1, 1, target_length, self.length);
        let mut forward_matrix = DpMatrixSparse::new(target_length, self.length, &row_bounds);

        // first we are going to generate n random
        // sequences drawn from the background
        // distribution and compute their forward scores
        let mut scores = vec![0.0; n];
        let mut rng = Pcg64::seed_from_u64(0);
        (0..n).for_each(|seq_idx| {
            forward_matrix.reuse(target_length, self.length, &row_bounds);
            let seq = Sequence::random_amino(target_length, &mut rng);

            // **NOTE: HMMER uses multi-hit mode to calibrate
            //         Tau, but we are using uni-hit mode
            let forward_score_nats = forward(self, &seq, &mut forward_matrix, &row_bounds);

            let null_score_nats = null_one_score(target_length);

            scores[seq_idx] = (forward_score_nats - null_score_nats).to_bits().value();
        });

        /// This is some black magic from a textbook:
        ///   Statistical Models and Methods for
        ///   Lifetime Data by Joseph F. Lawless
        ///   
        /// From HMMER:
        ///   Equation 4.1.6 from [Lawless82], pg. 143, and
        ///   its first derivative with respect to lambda,
        ///   for finding the ML fit to Gumbel lambda parameter.
        ///   This equation gives a result of zero for the maximum
        ///   likelihood lambda.
        fn lawless416(samples: &[f32], lambda: f32) -> (f32, f32) {
            // e_sum is the sum of e^(-lambda x_i)
            let mut e_sum = 0.0f32;
            // x_sum is the sum of x_i
            let mut x_sum = 0.0f32;
            // xe_sum is the sum of x_i * e^(-lambda x_i)
            let mut xe_sum = 0.0f32;
            // xe_sum is the sum of x_i^2 * e^(-lambda x_i)
            let mut xxe_sum = 0.0f32;

            samples.iter().for_each(|x| {
                e_sum += (-lambda * x).exp();
                x_sum += x;
                xe_sum += x * (-lambda * x).exp();
                xxe_sum += x.powi(2) * (-lambda * x).exp();
            });

            let fx = (1.0 / lambda) - (x_sum / samples.len() as f32) + (xe_sum / e_sum);
            let dfx = (xe_sum / e_sum).powi(2) - (xxe_sum / e_sum) - (1.0 / (lambda.powi(2)));

            (fx, dfx)
        }

        // now make an initial guess at lambda
        //
        // from hmmer:
        //   (Evans/Hastings/Peacock, Statistical Distributions, 2000, p.86)
        let sum: f32 = scores.iter().sum();
        let squared_sum: f32 = scores.iter().map(|s| s.powi(2)).sum();
        let sample_variance: f32 =
            (squared_sum - sum * sum / scores.len() as f32) / (scores.len() as f32 - 1.0);
        let mut gumbel_lambda: f32 = std::f32::consts::PI / (6.0 / sample_variance).sqrt();

        // now we do this Newton/Raphson root finding
        // thing until we have converged on lambda
        let tolerance: f32 = 1e-5;
        let mut newton_raphson_success = false;
        for _ in 0..=100 {
            let (fx, dfx) = lawless416(&scores, gumbel_lambda);
            if fx.abs() < tolerance {
                newton_raphson_success = true;
                break;
            }
            gumbel_lambda -= fx / dfx;

            if gumbel_lambda <= 0.0 {
                gumbel_lambda = 0.001;
            }
        }

        if !newton_raphson_success {
            panic!("newton/raphson failed");
        }

        // this is apparently substituting into equation 4.1.5
        // from Lawless[82] to solve for the mu parameter
        let e_sum: f32 = scores.iter().map(|s| (-gumbel_lambda * s).exp()).sum();
        let gumbel_mu = -(e_sum / scores.len() as f32).ln() / gumbel_lambda;

        // now that we've fit the gumbel lambda, we are going to find our tau

        /// Calculates the inverse CDF for a Gumbel distribution
        /// with parameters <mu> and <lambda>. That is, returns
        /// the quantile <x> at which the CDF is <p>.
        fn gumbel_inverse_cdf(p: f32, mu: f32, lambda: f32) -> f32 {
            mu - ((-p.ln()).ln() / lambda)
        }

        // basically, there is (maybe?) no good method for fitting an exponential,
        // so instead we have fit a Gumbel that we are going to use to pick our tau
        //
        // from hmmer:
        //   Explanation of the eqn below: first find the x at which the Gumbel tail
        //   mass is predicted to be equal to tailp. Then back up from that x
        //   by log(tailp)/lambda to set the origin of the exponential tail to 1.0
        //   instead of tailp.
        self.forward_tau = gumbel_inverse_cdf(1.0 - tail_probability, gumbel_mu, gumbel_lambda)
            + (tail_probability.ln() / self.forward_lambda);
    }

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
            consensus_sequence_bytes_utf8: vec![UTF8_SPACE],
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
                .consensus_sequence_bytes_utf8
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
            Trace::S_STATE | Trace::T_STATE => 0.0,
            Trace::N_STATE => match state_to {
                Trace::B_STATE => {
                    self.special_transition_score(Profile::SPECIAL_N_IDX, Profile::SPECIAL_MOVE_IDX)
                }
                Trace::N_STATE => {
                    self.special_transition_score(Profile::SPECIAL_N_IDX, Profile::SPECIAL_LOOP_IDX)
                }
                _ => panic!(),
            },
            Trace::B_STATE => match state_to {
                Trace::M_STATE => self.transition_score(Profile::BEGIN_TO_MATCH_IDX, idx_to - 1),
                _ => panic!(),
            },
            Trace::M_STATE => match state_to {
                Trace::M_STATE => self.transition_score(Profile::MATCH_TO_MATCH_IDX, idx_from),
                Trace::I_STATE => self.transition_score(Profile::MATCH_TO_INSERT_IDX, idx_from),
                Trace::D_STATE => self.transition_score(Profile::MATCH_TO_DELETE_IDX, idx_from),
                Trace::E_STATE => 0.0,
                _ => panic!(),
            },
            Trace::D_STATE => match state_to {
                Trace::M_STATE => self.transition_score(Profile::DELETE_TO_MATCH_IDX, idx_from),
                Trace::D_STATE => self.transition_score(Profile::DELETE_TO_DELETE_IDX, idx_from),
                Trace::E_STATE => 0.0,
                _ => panic!(),
            },
            Trace::I_STATE => match state_to {
                Trace::M_STATE => self.transition_score(Profile::INSERT_TO_MATCH_IDX, idx_from),
                Trace::I_STATE => self.transition_score(Profile::INSERT_TO_INSERT_IDX, idx_from),
                _ => panic!(),
            },
            Trace::E_STATE => match state_to {
                Trace::C_STATE => {
                    self.special_transition_score(Profile::SPECIAL_E_IDX, Profile::SPECIAL_MOVE_IDX)
                }
                Trace::J_STATE => {
                    self.special_transition_score(Profile::SPECIAL_E_IDX, Profile::SPECIAL_LOOP_IDX)
                }
                _ => panic!(),
            },
            Trace::J_STATE => match state_to {
                Trace::B_STATE => {
                    self.special_transition_score(Profile::SPECIAL_J_IDX, Profile::SPECIAL_MOVE_IDX)
                }
                Trace::J_STATE => {
                    self.special_transition_score(Profile::SPECIAL_J_IDX, Profile::SPECIAL_LOOP_IDX)
                }
                _ => panic!(),
            },
            Trace::C_STATE => match state_to {
                Trace::T_STATE => {
                    self.special_transition_score(Profile::SPECIAL_C_IDX, Profile::SPECIAL_MOVE_IDX)
                }
                Trace::C_STATE => {
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::structs::Sequence;

    #[test]
    fn test_calibrate_tau() -> anyhow::Result<()> {
        let mut seq = Sequence::from_utf8(
            concat!(
                "GNLLVILVILRNKKLRTPTNIFLLNLAVADLLVLLLVLPFSLVYALLEGDWVFGEVLCKL",
                "VTALDVVNLTASILLLTAISIDRYLAIVKPLKYKRIRTKRRALVLILVVWVLALLLSLPP",
                "LLFSGTKTESAEKEETVCLIDFPEEESTWEVSYTLLLSVLGFLLPLLVILVCYVRILRTL",
                "RKSAKKEKSRKKKSARKERKALKTLLVVVVVFVLCWLPYFILLLLDSLLKECESEKLVET",
                "ALLITLLLAYVNSCLNPIIY"
            )
            .as_bytes(),
        )?;
        seq.name = ">7tm_1-consensus".to_string();

        let hmm = Hmm::from_blosum_62_and_sequence(&seq)?;
        let mut profile = Profile::new(&hmm);

        profile.calibrate_tau(200, 100, 0.04);

        let correct_tau = -4.193134f32;
        let diff = (correct_tau - profile.forward_tau).abs();
        assert!(diff <= 1e-5);

        Ok(())
    }
}
