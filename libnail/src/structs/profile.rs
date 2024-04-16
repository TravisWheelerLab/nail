use crate::align::structs::{DpMatrixFlat, DpMatrixSparse, RowBounds, Trace};
use crate::align::{forward, length_bias_score};
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
use crate::util::{f32_vec_argmax, LogAbuse};

use std::fmt;
use std::fmt::Formatter;

use super::Sequence;

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

    pub fn calibrate_tau(&mut self, n: usize, target_length: usize) {
        self.configure_for_target_length(target_length);

        let mut row_bounds = RowBounds::new(target_length);
        row_bounds.fill_rectangle(1, 1, target_length, self.length);
        let mut forward_matrix = DpMatrixSparse::new(target_length, self.length, &row_bounds);

        let mut scores = vec![0.0; n];
        (0..n).for_each(|seq_idx| {
            forward_matrix.reuse(target_length, self.length, &row_bounds);
            let seq = Sequence::random_amino(target_length);
            let forward_score_nats = forward(self, &seq, &mut forward_matrix, &row_bounds);
            let null_score_nats = length_bias_score(target_length);

            let forward_score_bits =
                (forward_score_nats - null_score_nats) / std::f32::consts::LN_2;

            scores[seq_idx] = forward_score_bits;
        });

        //let scores = [
        //    -1.504366, 0.895846, -0.172907, -2.620571, -2.868223, -3.176270, -3.728831, -0.185057,
        //    -3.539029, -2.430213, -3.994046, -2.577063, -3.359447, -2.216906, -2.557495, -2.467121,
        //    -3.420379, -2.475419, -0.617509, -1.572560, -2.625687, -3.467364, -2.532328, -1.486470,
        //    -3.099069, -1.088266, -1.554731, -2.355674, -3.488132, -2.104149, -1.346401, -1.475866,
        //    -2.891736, -2.777906, -0.928439, -3.422777, -1.885561, -3.030625, -3.083956, -2.576361,
        //    -2.931466, -2.590394, -2.909705, -3.767259, -2.019206, -2.598923, -0.033133, -1.282793,
        //    1.145949, -2.516939, -3.042393, 0.450592, -1.811978, -3.199646, -2.705153, -2.248631,
        //    -3.127073, -2.809953, -2.711109, -2.604893, -1.082016, -3.801274, -3.265915, -2.958324,
        //    -3.936880, -2.592531, -1.580674, -3.264680, -2.049550, -2.643699, -2.725697, -3.362312,
        //    -3.266914, -3.090422, -3.677707, -2.412420, -2.727616, -2.725939, -3.437995, -3.088321,
        //    -3.474075, -0.222862, -3.244580, -2.778144, -3.807859, -2.530046, -3.687471, -3.098586,
        //    -0.062036, -2.713214, -1.529500, -3.595964, -3.413503, -2.127479, -2.375750, -3.464004,
        //    -3.203416, -1.626595, -3.734239, -2.022500, -2.140148, -1.103836, -2.039302, -2.788296,
        //    -2.467240, -3.288251, -1.950479, -2.874652, -3.324462, -3.470531, -1.925927, -3.411441,
        //    -3.593416, -2.867798, -1.993794, -3.205558, -2.991038, -2.961757, -1.953620, -3.768302,
        //    -2.976820, -1.868483, -3.590745, -3.227986, -3.009086, -2.223046, -2.981594, 1.806790,
        //    -3.989222, -0.917316, -2.295252, -3.476153, -0.518589, -2.940001, -3.148827, -2.876527,
        //    -2.515204, -2.406154, -2.913562, -2.836409, -2.924591, -2.777573, -3.988820, -3.669911,
        //    -3.195895, -3.341163, -2.802616, -1.503576, -1.748802, -3.406860, -2.448678, -2.638261,
        //    -3.303933, -3.765265, -2.579837, -2.350813, -1.209413, -3.172528, -3.000526, -3.602516,
        //    -2.859212, -3.444830, -3.297152, -2.862897, -2.368723, -1.845725, -1.768618, -2.923238,
        //    1.296241, -3.842624, -3.306856, -1.713498, -2.379782, -2.994811, -2.222899, 0.877470,
        //    -3.197263, -3.331733, -3.180032, -1.565291, 0.363498, -2.070105, -1.979063, -2.417173,
        //    -4.005939, -4.027726, -3.326597, -2.862119, -2.628895, -3.245254, -3.533534, -2.602711,
        //    -1.578742, -1.513021, -3.185244, -3.122952, -2.299559, -0.400301, -3.316555, -2.450500,
        //];

        fn lawless(samples: &[f32], lambda: f32) -> (f32, f32) {
            // e_sum is the sum of e^(-lambda x_i)
            let mut e_sum = 0.0;
            // x_sum is the sum of x_i
            let mut x_sum = 0.0;
            // xe_sum is the sum of x_i * e^(-lambda x_i)
            let mut xe_sum = 0.0;
            // xe_sum is the sum of x_i^2 * e^(-lambda x_i)
            let mut xxe_sum = 0.0;

            samples.iter().for_each(|x| {
                e_sum += (-lambda * x).exp();
                x_sum += x;
                xe_sum += x * (-lambda * x).exp();
                xxe_sum += x.powi(2) + (-lambda * x).exp();
            });

            let fx = (1.0 / lambda) - (x_sum / samples.len() as f32) + (xe_sum / e_sum);
            let dfx = (xe_sum / e_sum).powi(2) - (xxe_sum / e_sum) - (1.0 / (lambda.powi(2)));

            (fx, dfx)
        }

        let sum: f32 = scores.iter().sum();
        let squared_sum: f32 = scores.iter().map(|s| s.powi(2)).sum();

        let sample_variance: f32 =
            (squared_sum - sum * sum / scores.len() as f32) / (scores.len() as f32 - 1.0);
        let mut lambda: f32 = std::f32::consts::PI / (6.0 / sample_variance).sqrt();

        let tolerance: f32 = 1e-5;
        let mut newton_raphson_success = false;
        for _ in 0..=100 {
            let (fx, dfx) = lawless(&scores, lambda);
            println!("{fx}");

            if fx.abs() < tolerance {
                newton_raphson_success = true;
                break;
            }
            lambda -= fx / dfx;
            lambda = lambda.max(0.001);
            if lambda <= 0.0 {
                lambda = 0.001;
            }
        }
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

        profile.calibrate_tau(200, 100);
        panic!();
        Ok(())
    }
}
