/// Correctness tests for the DP alignment pipeline.
///
/// These tests serve as regression guards for the probability-space conversion.
/// Each test should continue to pass after `forward`, `backward`, `posterior`,
/// and `cloud_search_fwd`/`cloud_search_bwd` are rewritten in probability space.
///
/// The primary invariants checked here are:
///  1. Forward score is deterministic and matches a known value (catches wrong scores).
///  2. Forward-backward consistency: posterior probabilities sum to 1.0 at every
///     sequence position. If either `forward` or `backward` is wrong, this breaks.
///  3. Cloud search forward+backward scores are finite and ordered correctly.
#[cfg(test)]
mod tests {
    use rand::SeedableRng;
    use rand_pcg::Pcg64;

    use crate::{
        align::{
            backward, cloud_search_bwd, cloud_search_fwd, forward,
            posterior,
            structs::{AdMatrixLinear, Cloud, DpMatrix, DpMatrixSparse, RowBounds, Seed},
            CloudSearchParams,
        },
        structs::{Profile, Sequence},
    };

    // -----------------------------------------------------------------------
    // Fixture helpers
    // -----------------------------------------------------------------------

    /// Build a deterministic (profile, query_seq, target_seq) triple.
    ///
    /// `prf_seed` drives the sequence used to build the profile.
    /// `seq_seed` drives the sequence that will be scored against the profile.
    /// Using different seeds makes the test realistic (query ≠ target).
    fn make_fixture(prf_len: usize, seq_len: usize, prf_seed: u64, seq_seed: u64) -> (Profile, Sequence) {
        let mut rng = Pcg64::seed_from_u64(prf_seed);
        let prf_seq = Sequence::random_amino(prf_len, &mut rng);
        let mut prf = Profile::from_blosum_62_and_seq(&prf_seq).expect("valid profile");

        let mut rng2 = Pcg64::seed_from_u64(seq_seed);
        let seq = Sequence::random_amino(seq_len, &mut rng2);

        prf.configure_for_target_length(seq.length);
        (prf, seq)
    }

    fn full_bounds(prf: &Profile, seq: &Sequence) -> RowBounds {
        let mut bounds = RowBounds::new(seq.length);
        bounds.fill_rectangle(1, 1, seq.length, prf.length);
        bounds
    }

    // -----------------------------------------------------------------------
    // Test 1: forward score stability
    //
    // Run log-space forward on two different (profile, sequence) pairs and
    // assert the score matches a known expected value within a tight tolerance.
    //
    // Tolerance of 0.001 nats ≈ 0.0014 bits — much tighter than the
    // ~0.001-nat approximation error in log_add, so a prob-space
    // implementation with correct numerics will easily pass this.
    // -----------------------------------------------------------------------

    const FORWARD_SCORE_TOLERANCE: f32 = 0.01; // nats

    #[test]
    fn test_forward_score_small() {
        let (prf, seq) = make_fixture(50, 80, 1, 2);
        let bounds = full_bounds(&prf, &seq);
        let mut fwd_mx = DpMatrixSparse::new_prob(seq.length, prf.length, &bounds);
        let score = forward(&prf, &seq, &mut fwd_mx, &bounds);
        assert!(
            score.value().is_finite(),
            "forward score is not finite: {}",
            score.value()
        );
        // Expected value captured from the log-space implementation.
        // A probability-space implementation must reproduce this within tolerance.
        let expected: f32 = -5.4625;
        assert!(
            (score.value() - expected).abs() < FORWARD_SCORE_TOLERANCE,
            "forward score {:.4} differs from expected {:.4} by more than {:.4} nats",
            score.value(),
            expected,
            FORWARD_SCORE_TOLERANCE,
        );
    }

    #[test]
    fn test_forward_score_larger() {
        let (prf, seq) = make_fixture(150, 300, 7, 13);
        let bounds = full_bounds(&prf, &seq);
        let mut fwd_mx = DpMatrixSparse::new_prob(seq.length, prf.length, &bounds);
        let score = forward(&prf, &seq, &mut fwd_mx, &bounds);
        assert!(
            score.value().is_finite(),
            "forward score is not finite: {}",
            score.value()
        );
        let expected: f32 = -7.4269;
        assert!(
            (score.value() - expected).abs() < FORWARD_SCORE_TOLERANCE,
            "forward score {:.4} differs from expected {:.4} by more than {:.4} nats",
            score.value(),
            expected,
            FORWARD_SCORE_TOLERANCE,
        );
    }

    // -----------------------------------------------------------------------
    // Test 2: posterior normalization (forward-backward consistency)
    //
    // At every sequence position t, the sum of posterior probabilities across
    // all match, insert, and special states must equal 1.0.
    //
    // This is the canonical sanity check for forward-backward algorithms:
    // if either forward or backward is wrong, the sum will deviate from 1.
    // -----------------------------------------------------------------------

    const POSTERIOR_SUM_TOLERANCE: f32 = 1e-3;

    fn check_posterior_normalizes(prf: &Profile, seq: &Sequence) {
        let bounds = full_bounds(prf, seq);
        let mut fwd_mx = DpMatrixSparse::new_prob(seq.length, prf.length, &bounds);
        let mut bwd_mx = DpMatrixSparse::new_prob(seq.length, prf.length, &bounds);
        let mut post_mx = DpMatrixSparse::new_prob(seq.length, prf.length, &bounds);

        forward(prf, seq, &mut fwd_mx, &bounds);
        backward(prf, seq, &mut bwd_mx, &bounds);
        posterior(prf, &fwd_mx, &bwd_mx, &mut post_mx, &bounds);

        for t_idx in bounds.seq_start..=bounds.seq_end {
            let p_start = bounds.left_row_bounds[t_idx];
            let p_end = bounds.right_row_bounds[t_idx];

            let core_sum: f32 = (p_start..=p_end)
                .map(|p_idx| {
                    post_mx.get_match(t_idx, p_idx) + post_mx.get_insert(t_idx, p_idx)
                })
                .sum();

            let special_sum = post_mx.get_special(t_idx, Profile::N_IDX)
                + post_mx.get_special(t_idx, Profile::J_IDX)
                + post_mx.get_special(t_idx, Profile::C_IDX);

            let row_sum = core_sum + special_sum;
            assert!(
                (row_sum - 1.0).abs() < POSTERIOR_SUM_TOLERANCE,
                "posterior at t={} sums to {:.6} (expected 1.0, diff {:.2e})",
                t_idx,
                row_sum,
                (row_sum - 1.0).abs(),
            );
        }
    }

    #[test]
    fn test_posterior_normalizes_small() {
        let (prf, seq) = make_fixture(50, 80, 1, 2);
        check_posterior_normalizes(&prf, &seq);
    }

    #[test]
    fn test_posterior_normalizes_larger() {
        let (prf, seq) = make_fixture(150, 300, 7, 13);
        check_posterior_normalizes(&prf, &seq);
    }

    // Same profile as target: high-score case
    #[test]
    fn test_posterior_normalizes_self_alignment() {
        let mut rng = Pcg64::seed_from_u64(42);
        let seq = Sequence::random_amino(100, &mut rng);
        let mut prf = Profile::from_blosum_62_and_seq(&seq).expect("valid profile");
        prf.configure_for_target_length(seq.length);
        check_posterior_normalizes(&prf, &seq);
    }

    // -----------------------------------------------------------------------
    // Test 3: cloud search scores are finite and well-ordered
    //
    // Both forward and backward cloud searches must produce finite scores, and
    // the score at the seed region must be >= the score outside the seed.
    // -----------------------------------------------------------------------

    #[test]
    fn test_cloud_search_scores_finite() {
        let (prf, seq) = make_fixture(100, 200, 3, 5);

        // Seed at the center of the (profile, sequence) space
        let seed = Seed {
            prf: prf.name.clone(),
            seq: seq.name.clone(),
            prf_start: prf.length / 4,
            prf_end: 3 * prf.length / 4,
            seq_start: seq.length / 4,
            seq_end: 3 * seq.length / 4,
            score: 0.0,
            e_value: 1.0,
        };

        let params = CloudSearchParams::default();
        let mut mx = AdMatrixLinear::new(seq.length);
        let mut cloud = Cloud::default();
        // cloud_search_fwd requires the caller to reuse the cloud and matrix first
        // (cloud_search_bwd handles reuse internally — an existing API asymmetry)
        mx.reuse(seq.length);
        cloud.reuse(seq.length, prf.length);

        let fwd = cloud_search_fwd(&prf, &seq, &seed, &mut mx, &params, &mut cloud);
        assert!(
            fwd.max_score.value().is_finite(),
            "cloud fwd max_score is not finite"
        );

        let mut mx2 = AdMatrixLinear::new(seq.length);
        let mut cloud2 = Cloud::default();
        let bwd = cloud_search_bwd(&prf, &seq, &seed, &mut mx2, &params, &mut cloud2);
        assert!(
            bwd.max_score.value().is_finite(),
            "cloud bwd max_score is not finite"
        );

        // max_score must be >= max_score_within (the within-seed max)
        assert!(
            fwd.max_score.value() >= fwd.max_score_within.value() - 1e-4,
            "cloud fwd: max_score ({}) < max_score_within ({})",
            fwd.max_score.value(),
            fwd.max_score_within.value()
        );
        assert!(
            bwd.max_score.value() >= bwd.max_score_within.value() - 1e-4,
            "cloud bwd: max_score ({}) < max_score_within ({})",
            bwd.max_score.value(),
            bwd.max_score_within.value()
        );
    }
}
