use crate::align::structs::{DpMatrix, RowBounds};
use crate::log_sum;
use crate::structs::{Profile, Sequence};
use crate::util::log_add;

use super::CloudSearchResults;

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
    fn min(self, other: Self) -> Self;
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

    fn min(self, other: Self) -> Self {
        Nats(self.0.min(other.0))
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

    fn min(self, other: Self) -> Self {
        Bits(self.0.min(other.0))
    }
}

pub fn p_value(score: impl Score, lambda: f32, tau: f32) -> f64 {
    (-lambda as f64 * ((score.bits().value()) as f64 - tau as f64)).exp()
}

pub fn e_value(p_value: f64, num_targets: usize) -> f64 {
    p_value * num_targets as f64
}

/// Compute the cloud score: the approximation of the forward score of the entire cloud.
pub fn cloud_score(
    forward_scores: &CloudSearchResults,
    reverse_scores: &CloudSearchResults,
) -> Nats {
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
    // TODO: prevent these allocations?
    let mut expected_prob_ratios: Vec<f32> = vec![0.0; Profile::MAX_DEGENERATE_ALPHABET_SIZE];
    let mut match_sums: Vec<f32> = vec![0.0; profile.length + 1];
    let mut insert_sums: Vec<f32> = vec![0.0; profile.length + 1];
    let mut core_posteriors: Vec<f32> = vec![0.0; target.length + 1];
    let mut core_state_sum: f32 = 0.0;

    // what: for each position in the model, take the sum of
    //       the posteriors in the match and insert state
    //
    // why: this gives us the expected number of times that
    //      each state was used in generating the target sequence
    //
    //      for example:
    //
    //              P_1   P_2   P_3    N     C
    //
    //    T_1   M   0.05  0.05  0.00  0.80  0.00
    //          I   0.10  0.00  0.00
    //          D   0.00  0.00  0.00
    //
    //    T_2   M   0.30  0.20  0.10  0.40  0.00
    //          I   0.00  0.00  0.00
    //          D   0.00  0.00  0.00
    //
    //    T_3   M   0.25  0.50  0.25  0.00  0.00
    //          I   0.00  0.00  0.00
    //          D   0.00  0.00  0.00
    //
    //    T_4   M   0.10  0.20  0.30  0.00  0.40
    //          I   0.00  0.00  0.00
    //          D   0.00  0.00  0.00
    //
    //    T_5   M   0.00  0.05  0.05  0.00  0.80
    //          I   0.00  0.00  0.10
    //          D   0.00  0.00  0.00
    //
    //        ----------------------------------
    //          M   0.70  1.00  0.70  1.20  1.20
    //          I   0.10  0.00  0.10
    //
    //                ^     ^     ^     ^     ^
    //              these are the expected numbers
    //              of times each state is used
    //
    for target_idx in row_bounds.target_start..=row_bounds.target_end {
        let profile_start_in_current_row = row_bounds.left_row_bounds[target_idx];
        let profile_end_in_current_row = row_bounds.right_row_bounds[target_idx];

        for profile_idx in profile_start_in_current_row..profile_end_in_current_row {
            match_sums[profile_idx] += posterior_matrix.get_match(target_idx, profile_idx);
            insert_sums[profile_idx] += posterior_matrix.get_insert(target_idx, profile_idx);
        }

        // the posterior probability of being in a special
        // state at this target position is the sum of
        // the individual posteriors of each special state
        let special_posterior = posterior_matrix.get_special(target_idx, Profile::SPECIAL_N_IDX)
            + posterior_matrix.get_special(target_idx, Profile::SPECIAL_J_IDX)
            + posterior_matrix.get_special(target_idx, Profile::SPECIAL_C_IDX);

        let core_posterior = 1.0 - special_posterior;
        core_posteriors[target_idx] = core_posterior;
        core_state_sum += core_posterior;
    }

    // now that we have the expected number of state usages,
    // we are going to compute the expected probability
    // ratios that we use to determine our composition
    // bias score adjustment (i.e. null two score)
    expected_prob_ratios
        .iter_mut()
        .enumerate()
        .take(Profile::MAX_ALPHABET_SIZE)
        .for_each(|(residue, ratio)| {
            for profile_idx in 1..profile.length {
                let match_contribution =
                    match_sums[profile_idx] * profile.match_score(residue, profile_idx).exp();

                let insert_contribution =
                    insert_sums[profile_idx] * profile.insert_score(residue, profile_idx).exp();

                *ratio += match_contribution + insert_contribution;
            }
            let match_contribution =
                match_sums[profile.length] * profile.match_score(residue, profile.length).exp();

            *ratio += match_contribution;
            *ratio /= core_state_sum;
        });

    // we set the scores for the degenerate characters to the
    // average of the scores of residues that they may represent
    // for example:
    //   the degenerate character B may either be a D or an N, so
    //   the score for B is the average of the scores for D and N

    // B ->  [D, N]
    // 21 -> [2, 11]
    expected_prob_ratios[21] = (expected_prob_ratios[2] + expected_prob_ratios[11]) / 2.0;
    // J ->  [I, L]
    // 22 -> [7, 9]
    expected_prob_ratios[22] = (expected_prob_ratios[7] + expected_prob_ratios[9]) / 2.0;
    // Z ->  [E, Q]
    // 23 -> [3, 13]
    expected_prob_ratios[23] = (expected_prob_ratios[3] + expected_prob_ratios[13]) / 2.0;
    // U ->  [C]
    // 24 -> [1]
    expected_prob_ratios[24] = expected_prob_ratios[1];
    // O ->  [K]
    // 25 -> [8]
    expected_prob_ratios[25] = expected_prob_ratios[8];
    // X ->  [any]
    // 26 -> [0..19]
    expected_prob_ratios[26] = expected_prob_ratios[0..20].iter().sum::<f32>() / 20.0;

    expected_prob_ratios[Profile::GAP_INDEX] = 1.0;
    expected_prob_ratios[Profile::NON_RESIDUE_IDX] = 1.0;
    expected_prob_ratios[Profile::MISSING_DATA_IDX] = 1.0;

    let expected_scores: Vec<_> = expected_prob_ratios.into_iter().map(|r| r.ln()).collect();

    let mut null_two_score = 0.0;
    (row_bounds.target_start..=row_bounds.target_end)
        .map(|idx| {
            (
                core_posteriors[idx],
                expected_scores[target.digital_bytes[idx] as usize],
            )
        })
        .for_each(|(core_posterior, expected_score)| {
            // we weight each residue's contribution to the
            // bias by the posterior probability that
            // the residue was emitted by a core model state
            null_two_score += core_posterior * expected_score;
        });

    // this is "omega" in hmmer
    //   essentially, we have a strong prior expecation that
    //   the sequence does not have a composition bias
    let null_two_prior = (1.0f32 / 256.0).ln();
    null_two_score = log_sum!(0.0, null_two_prior + null_two_score);
    Nats(null_two_score)
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
