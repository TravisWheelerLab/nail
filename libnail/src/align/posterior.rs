use crate::align::structs::{DpMatrix, DpMatrixSparse, RowBounds};
use crate::structs::Profile;

/// Probability-space posterior decoding.
///
/// Computes posterior probabilities from probability-space forward and backward matrices.
/// For core states: posterior[t][p] ∝ fwd[t][p] * bwd[t][p]
/// For special states (N, C, J): posterior[t] ∝ fwd[t-1] * loop * bwd[t]
///   These need a correction factor (forward row scale at t) because they
///   use fwd[t-1] which has a different cumulative scale than fwd[t].
///
/// All posteriors are normalized so that each row sums to 1.0.
pub fn posterior(
    profile: &Profile,
    forward_matrix: &DpMatrixSparse,
    backward_matrix: &DpMatrixSparse,
    posterior_matrix: &mut DpMatrixSparse,
    row_bounds: &RowBounds,
) {
    posterior_matrix.set_special(row_bounds.seq_start - 1, Profile::E_IDX, 0.0);
    posterior_matrix.set_special(row_bounds.seq_start - 1, Profile::N_IDX, 0.0);
    posterior_matrix.set_special(row_bounds.seq_start - 1, Profile::J_IDX, 0.0);
    posterior_matrix.set_special(row_bounds.seq_start - 1, Profile::B_IDX, 0.0);
    posterior_matrix.set_special(row_bounds.seq_start - 1, Profile::C_IDX, 0.0);

    let profile_start_in_first_row = row_bounds.left_row_bounds[row_bounds.seq_start];
    let profile_end_in_first_row = row_bounds.right_row_bounds[row_bounds.seq_start];

    for profile_idx in (profile_start_in_first_row - 1)..=profile_end_in_first_row {
        posterior_matrix.set_match(row_bounds.seq_start - 1, profile_idx, 0.0);
        posterior_matrix.set_insert(row_bounds.seq_start - 1, profile_idx, 0.0);
        posterior_matrix.set_delete(row_bounds.seq_start - 1, profile_idx, 0.0);
    }

    for target_idx in row_bounds.seq_start..=row_bounds.seq_end {
        let lp = row_bounds.left_row_bounds[target_idx];
        let rp = row_bounds.right_row_bounds[target_idx];

        posterior_matrix.set_match(target_idx, lp - 1, 0.0);
        posterior_matrix.set_insert(target_idx, lp - 1, 0.0);
        posterior_matrix.set_delete(target_idx, lp - 1, 0.0);

        // Element-wise product for M, I, D over [lp, rp]
        {
            let fwd_m = forward_matrix.m_slice(target_idx, lp, rp);
            let bwd_m = backward_matrix.m_slice(target_idx, lp, rp);
            let fwd_i = forward_matrix.i_slice(target_idx, lp, rp);
            let bwd_i = backward_matrix.i_slice(target_idx, lp, rp);
            let (post_m, post_i, post_d) = posterior_matrix.core_slices_mut(target_idx, lp, rp);
            post_m
                .iter_mut()
                .zip(fwd_m.iter().zip(bwd_m.iter()))
                .for_each(|(pm, (fm, bm))| *pm = fm * bm);
            post_i
                .iter_mut()
                .zip(fwd_i.iter().zip(bwd_i.iter()))
                .for_each(|(pi, (fi, bi))| *pi = fi * bi);
            post_d.fill(0.0);
        }

        // Zero insert at last profile position if at end of profile
        if rp == profile.length {
            posterior_matrix.set_insert(target_idx, rp, 0.0);
        }

        // Denominator from core states
        let mut denominator: f32 = {
            let (post_m, post_i, _) = posterior_matrix.core_slices_mut(target_idx, lp, rp);
            post_m.iter().sum::<f32>() + post_i.iter().sum::<f32>()
        };

        posterior_matrix.set_special(target_idx, Profile::E_IDX, 0.0);

        // Special states: these use fwd[t-1] which has cumulative scale S^fwd_{t-1},
        // while core states use fwd[t] with scale S^fwd_t. The ratio is s^fwd_t
        // (the per-row scale factor for forward row t). We must multiply the special
        // state posteriors by 1/s^fwd_t to put them on the same scale as core states,
        // OR equivalently multiply core states by s^fwd_t. Since we normalize anyway,
        // we just need them all on the same relative scale. We multiply specials by
        // 1/fwd_scale_t to compensate.
        let fwd_scale = forward_matrix.get_row_scale(target_idx);
        let scale_correction = if fwd_scale > 0.0 { 1.0 / fwd_scale } else { 0.0 };

        // N state: P(emit t via N) ∝ fwd_N[t-1] * loop_N * bwd_N[t]
        let n_post = forward_matrix.get_special(target_idx - 1, Profile::N_IDX)
            * profile.special_transition_prob(Profile::N_IDX, Profile::SPECIAL_LOOP_IDX)
            * backward_matrix.get_special(target_idx, Profile::N_IDX)
            * scale_correction;
        posterior_matrix.set_special(target_idx, Profile::N_IDX, n_post);

        // J state: P(emit t via J) ∝ fwd_J[t-1] * loop_J * bwd_J[t]
        let j_post = forward_matrix.get_special(target_idx - 1, Profile::J_IDX)
            * profile.special_transition_prob(Profile::J_IDX, Profile::SPECIAL_LOOP_IDX)
            * backward_matrix.get_special(target_idx, Profile::J_IDX)
            * scale_correction;
        posterior_matrix.set_special(target_idx, Profile::J_IDX, j_post);

        posterior_matrix.set_special(target_idx, Profile::B_IDX, 0.0);

        // C state: P(emit t via C) ∝ fwd_C[t-1] * loop_C * bwd_C[t]
        let c_post = forward_matrix.get_special(target_idx - 1, Profile::C_IDX)
            * profile.special_transition_prob(Profile::C_IDX, Profile::SPECIAL_LOOP_IDX)
            * backward_matrix.get_special(target_idx, Profile::C_IDX)
            * scale_correction;
        posterior_matrix.set_special(target_idx, Profile::C_IDX, c_post);

        denominator += n_post + j_post + c_post;

        // Normalize (guard against zero denominator from empty rows)
        if denominator <= 0.0 {
            continue;
        }
        let inv = 1.0 / denominator;

        // Normalize M over [lp, rp]; normalize I over [lp, rp-1].
        // I at rp is excluded from normalization (matching original behavior:
        // it was either 0.0 or contributed to denominator but not re-scaled).
        {
            let (post_m, post_i, _) = posterior_matrix.core_slices_mut(target_idx, lp, rp);
            post_m.iter_mut().for_each(|v| *v *= inv);
            let i_main_len = post_i.len() - 1;
            post_i[..i_main_len].iter_mut().for_each(|v| *v *= inv);
        }

        posterior_matrix.set_special(target_idx, Profile::N_IDX, n_post * inv);
        posterior_matrix.set_special(target_idx, Profile::J_IDX, j_post * inv);
        posterior_matrix.set_special(target_idx, Profile::C_IDX, c_post * inv);
    }
}
