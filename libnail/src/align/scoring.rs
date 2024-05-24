use crate::align::structs::{DpMatrix, RowBounds};
use crate::log_sum;
use crate::structs::{Profile, Sequence};
use crate::util::{log_add, LogAbuse};

use super::CloudSearchScores;

/// A wrapper around f32 to describe nats
#[derive(Clone, Copy)]
pub struct Nats(pub f32);
impl Nats {
    pub fn value(&self) -> f32 {
        self.0
    }

    pub fn to_bits(self) -> Bits {
        Bits(self.0 / std::f32::consts::LN_2)
    }
}

impl std::fmt::Debug for Nats {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "Nats({})", self.0)
    }
}

impl std::ops::Add for Nats {
    type Output = Nats;

    fn add(self, rhs: Self) -> Self::Output {
        Nats(self.0 + rhs.0)
    }
}

impl std::ops::Sub for Nats {
    type Output = Nats;

    fn sub(self, rhs: Self) -> Self::Output {
        Nats(self.0 - rhs.0)
    }
}

impl std::ops::Add<Bits> for Nats {
    type Output = Nats;

    fn add(self, rhs: Bits) -> Self::Output {
        self + rhs.to_nats()
    }
}

impl std::ops::Sub<Bits> for Nats {
    type Output = Nats;

    fn sub(self, rhs: Bits) -> Self::Output {
        self - rhs.to_nats()
    }
}

/// A wrapper around f32 to describe bits
#[derive(Clone, Copy)]
pub struct Bits(pub f32);
impl Bits {
    pub fn value(&self) -> f32 {
        self.0
    }

    pub fn to_nats(self) -> Nats {
        Nats(self.0 * std::f32::consts::LN_2)
    }
}

impl std::fmt::Debug for Bits {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "Bits({})", self.0)
    }
}

impl std::ops::Add for Bits {
    type Output = Bits;

    fn add(self, rhs: Self) -> Self::Output {
        Bits(self.0 + rhs.0)
    }
}

impl std::ops::Sub for Bits {
    type Output = Bits;

    fn sub(self, rhs: Self) -> Self::Output {
        Bits(self.0 - rhs.0)
    }
}

impl std::ops::Add<Nats> for Bits {
    type Output = Bits;

    fn add(self, rhs: Nats) -> Self::Output {
        self + rhs.to_bits()
    }
}

impl std::ops::Sub<Nats> for Bits {
    type Output = Bits;

    fn sub(self, rhs: Nats) -> Self::Output {
        self - rhs.to_bits()
    }
}

/// A trait to generically accept a score value as either bits or nats
pub trait Score {
    fn bits(self) -> Bits;
    fn nats(self) -> Nats;
    fn max(self, other: Self) -> Self;
}

impl Score for Nats {
    fn bits(self) -> Bits {
        self.to_bits()
    }

    fn nats(self) -> Nats {
        self
    }

    fn max(self, other: Self) -> Self {
        Nats(self.0.max(other.0))
    }
}

impl Score for Bits {
    fn bits(self) -> Bits {
        self
    }

    fn nats(self) -> Nats {
        self.to_nats()
    }

    fn max(self, other: Self) -> Self {
        Bits(self.0.max(other.0))
    }
}

pub fn p_value(score: impl Score, lambda: f32, tau: f32) -> f64 {
    (-lambda as f64 * ((score.bits().value()) as f64 - tau as f64)).exp()
}

pub fn e_value(p_value: f64, num_targets: usize) -> f64 {
    p_value * num_targets as f64
}

/// Compute the cloud score: the approximation of the forward score of the entire cloud.
pub fn cloud_score(forward_scores: &CloudSearchScores, reverse_scores: &CloudSearchScores) -> Nats {
    // this approximates the score for the forward
    // cloud that extends past the seed end point
    let disjoint_forward_score = forward_scores.max_score - forward_scores.max_score_within;

    // this approximates the score for the reverse
    // cloud that extends past the seed start point
    let disjoint_reverse_score = reverse_scores.max_score - reverse_scores.max_score_within;

    // this approximates the score of the intersection
    // of the forward and reverse clouds
    let intersection_score = forward_scores
        .max_score_within
        .max(reverse_scores.max_score_within);

    intersection_score + disjoint_forward_score + disjoint_reverse_score
}

/// Compute the null one score adjustment: the sum of the background
/// transitions across the length of the target sequence.
pub fn null_one_score(target_length: usize) -> Nats {
    let p1 = (target_length as f32) / (target_length as f32 + 1.0);
    Nats(target_length as f32 * p1.ln() + (1.0 - p1).ln())
}

/// Compute the null two score adjustment: the composition bias.
pub fn null_two_score(
    posterior_matrix: &impl DpMatrix,
    profile: &Profile,
    target: &Sequence,
    row_bounds: &RowBounds,
) -> Nats {
    let sub_target_length = (row_bounds.target_end - row_bounds.target_start + 1) as f32;

    // TODO: need to prevent these allocations
    //       the strategy used in hmmer (using the first row in the PP-DP matrix)
    //       won't work here, since we have a sparse matrix
    let mut match_values: Vec<f32> = vec![0.0; profile.length + 1];
    let mut insert_values: Vec<f32> = vec![0.0; profile.length + 1];

    let mut n_value: f32 = 0.0;
    let mut j_value: f32 = 0.0;
    let mut c_value: f32 = 0.0;

    // calculate the expected number of times that each emitting
    // state was used in generating the residues in this domain
    for target_idx in row_bounds.target_start..=row_bounds.target_end {
        let profile_start_in_current_row = row_bounds.left_row_bounds[target_idx];
        let profile_end_in_current_row = row_bounds.right_row_bounds[target_idx];

        for profile_idx in profile_start_in_current_row..profile_end_in_current_row {
            // for profile_idx in 1..=profile.length {
            match_values[profile_idx] += posterior_matrix.get_match(target_idx, profile_idx);
            insert_values[profile_idx] += posterior_matrix.get_insert(target_idx, profile_idx);
        }

        n_value += posterior_matrix.get_special(target_idx, Profile::SPECIAL_N_IDX);
        j_value += posterior_matrix.get_special(target_idx, Profile::SPECIAL_J_IDX);
        c_value += posterior_matrix.get_special(target_idx, Profile::SPECIAL_C_IDX);
    }

    // convert the expected numbers to log frequencies where:
    //    the numerator is the expected number
    //    the denominator is the length of the part of the target that was aligned
    match_values
        .iter_mut()
        .for_each(|v| *v = v.ln_or_max() - sub_target_length.ln());

    insert_values
        .iter_mut()
        .for_each(|v| *v = v.ln_or_max() - sub_target_length.ln());

    n_value = n_value.ln_or_max() - sub_target_length.ln();
    j_value = j_value.ln_or_max() - sub_target_length.ln();
    c_value = c_value.ln_or_max() - sub_target_length.ln();

    // from hmmer:
    //   Calculate null2's log odds emission probabilities, by taking
    //   posterior weighted sum over all emission vectors used in paths
    //   explaining the domain.
    let mut null2: Vec<f32> = vec![-f32::INFINITY; Profile::MAX_DEGENERATE_ALPHABET_SIZE];
    let x_factor = log_sum!(n_value, j_value, c_value);

    for alphabet_idx in 0..Profile::MAX_ALPHABET_SIZE {
        for profile_idx in 1..profile.length {
            null2[alphabet_idx] = log_sum!(
                null2[alphabet_idx],
                match_values[profile_idx] + profile.match_score(alphabet_idx, profile_idx),
                insert_values[profile_idx] + profile.insert_score(alphabet_idx, profile_idx)
            );
        }
        null2[alphabet_idx] = log_sum!(
            null2[alphabet_idx],
            match_values[profile.length] + profile.match_score(alphabet_idx, profile.length),
            x_factor
        );
    }

    null2[0..20].iter_mut().for_each(|v| *v = v.exp());

    // from hmmer:
    //   now null2[x] = \frac{f_d(x)}{f_0(x)} for all x in alphabet,
    //   0..K-1, where f_d(x) are the ad hoc "null2" residue frequencies
    //   for this envelope.

    // we set the scores for the degenerate characters to the
    // average of the scores of residues that they may represent
    // for example:
    //   the degenerate character B may either be a D or an N, so
    //   the score for B is the average of the scores for D and N

    // B ->  [D, N]
    // 21 -> [2, 11]
    null2[21] = (null2[2] + null2[11]) / 2.0;
    // J ->  [I, L]
    // 22 -> [7, 9]
    null2[22] = (null2[7] + null2[9]) / 2.0;
    // Z ->  [E, Q]
    // 23 -> [3, 13]
    null2[23] = (null2[3] + null2[13]) / 2.0;
    // U ->  [C]
    // 24 -> [1]
    null2[24] = null2[1];
    // O ->  [K]
    // 25 -> [8]
    null2[25] = null2[8];
    // X ->  [any]
    // 26 -> [0..19]
    null2[26] = null2[0..20].iter().sum::<f32>() / 20.0;

    null2[Profile::GAP_INDEX] = 1.0;
    null2[Profile::NON_RESIDUE_IDX] = 1.0;
    null2[Profile::MISSING_DATA_IDX] = 1.0;

    let mut null2_score = 0.0;
    for residue in &target.digital_bytes[row_bounds.target_start..=row_bounds.target_end] {
        null2_score += null2[*residue as usize].ln();
    }

    // this is "omega" in hmmer
    // TODO: figure out more about this
    // TODO: also, do we really want to do this omega business here?
    let null2_prior = 1.0 / 256.0;
    null2_score = log_sum!(0.0, null2_prior + null2_score);

    Nats(null2_score)
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_nats_ops() -> anyhow::Result<()> {
        let a = Nats(10.0);
        let b = Nats(10.0);
        let c = a + b;
        assert_eq!(c.value(), 20.0);

        let a = Nats(20.0);
        let b = Nats(10.0);
        let c = a - b;
        assert_eq!(c.value(), 10.0);

        Ok(())
    }

    #[test]
    fn test_bits_ops() -> anyhow::Result<()> {
        let a = Bits(10.0);
        let b = Bits(10.0);
        let c = a + b;
        assert_eq!(c.value(), 20.0);

        let a = Bits(20.0);
        let b = Bits(10.0);
        let c = a - b;
        assert_eq!(c.value(), 10.0);

        Ok(())
    }

    #[test]
    fn test_nats_bits_conversion() -> anyhow::Result<()> {
        let tolerance = 1e-9;
        let nats = Nats(10.0f32.ln());
        let bits = Bits(10.0f32.log2());

        assert!((nats.value() - bits.to_nats().value()).abs() < tolerance);
        assert!((bits.value() - nats.to_bits().value()).abs() < tolerance);
        Ok(())
    }

    #[test]
    fn test_nats_bits_ops() -> anyhow::Result<()> {
        let a = Nats(10.0f32.ln());
        let b = Bits(10.0f32.log2());
        let c = a + b;
        let correct = a + a;
        assert_eq!(c.value(), correct.value());

        let a = Bits(10.0f32.log2());
        let b = Nats(10.0f32.ln());
        let c = a + b;
        let correct = a + a;
        assert_eq!(c.value(), correct.value());

        Ok(())
    }
}
