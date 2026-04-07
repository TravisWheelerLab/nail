use crate::align::structs::{DpMatrix, DpMatrixSparse, RowBounds};
use crate::structs::profile::Transition;
use crate::structs::{Profile, Sequence};

use super::Nats;

/// Probability-space Forward algorithm with per-row scaling.
///
/// All DP values are stored as probabilities (not log-probabilities).
/// Per-row scaling prevents underflow: after computing each row, all values
/// are divided by the row's maximum value, and log(max) is accumulated in `xsc`.
///
/// The final score in nats is: ln(C_stored[T]) + xsc + background_correction + C->T.
pub fn forward(
    profile: &Profile,
    target: &Sequence,
    dp_matrix: &mut DpMatrixSparse,
    bounds: &RowBounds,
) -> Nats {
    // Initial conditions in probability space
    // N[seq_start - 1] = 1.0 (we start in the N state with probability 1)
    dp_matrix.set_special(bounds.seq_start - 1, Profile::N_IDX, 1.0);
    dp_matrix.set_special(
        bounds.seq_start - 1,
        Profile::B_IDX,
        profile.special_transition_prob(Profile::N_IDX, Profile::SPECIAL_MOVE_IDX),
    );
    dp_matrix.set_special(bounds.seq_start - 1, Profile::E_IDX, 0.0);
    dp_matrix.set_special(bounds.seq_start - 1, Profile::C_IDX, 0.0);

    let mut xsc = 0.0f32; // accumulated log-scale

    for target_idx in bounds.seq_start..=bounds.seq_end {
        let residue = target.digital_bytes[target_idx] as usize;
        let b_prev = dp_matrix.get_special(target_idx - 1, Profile::B_IDX);

        // core_max tracked inline to avoid a separate read-only pass over the row.
        let mut core_max = 0.0f32;

        for profile_idx in bounds.left_row_bounds[target_idx]..bounds.right_row_bounds[target_idx]
        {
            // match state
            let m_val = (dp_matrix.get_match(target_idx - 1, profile_idx - 1)
                * profile.transition_prob(Transition::MM as usize, profile_idx - 1)
                + dp_matrix.get_insert(target_idx - 1, profile_idx - 1)
                    * profile.transition_prob(Transition::IM as usize, profile_idx - 1)
                + b_prev * profile.transition_prob(Transition::BM as usize, profile_idx - 1)
                + dp_matrix.get_delete(target_idx - 1, profile_idx - 1)
                    * profile.transition_prob(Transition::DM as usize, profile_idx - 1))
                * profile.match_prob(residue, profile_idx);
            dp_matrix.set_match(target_idx, profile_idx, m_val);

            // insert state
            let i_val = (dp_matrix.get_match(target_idx - 1, profile_idx)
                * profile.transition_prob(Transition::MI as usize, profile_idx)
                + dp_matrix.get_insert(target_idx - 1, profile_idx)
                    * profile.transition_prob(Transition::II as usize, profile_idx))
                * profile.insert_prob(residue, profile_idx);
            dp_matrix.set_insert(target_idx, profile_idx, i_val);

            // delete state
            let d_val = dp_matrix.get_match(target_idx, profile_idx - 1)
                * profile.transition_prob(Transition::MD as usize, profile_idx - 1)
                + dp_matrix.get_delete(target_idx, profile_idx - 1)
                    * profile.transition_prob(Transition::DD as usize, profile_idx - 1);
            dp_matrix.set_delete(target_idx, profile_idx, d_val);

            core_max = core_max.max(m_val).max(i_val);

            // E state (accumulate)
            let e_prev = dp_matrix.get_special(target_idx, Profile::E_IDX);
            dp_matrix.set_special(target_idx, Profile::E_IDX, e_prev + m_val + d_val);
        }

        // --- unrolled last column: profile_idx = right_row_bounds[target_idx] ---
        let profile_idx = bounds.right_row_bounds[target_idx];

        let m_val = (dp_matrix.get_match(target_idx - 1, profile_idx - 1)
            * profile.transition_prob(Transition::MM as usize, profile_idx - 1)
            + dp_matrix.get_insert(target_idx - 1, profile_idx - 1)
                * profile.transition_prob(Transition::IM as usize, profile_idx - 1)
            + b_prev * profile.transition_prob(Transition::BM as usize, profile_idx - 1)
            + dp_matrix.get_delete(target_idx - 1, profile_idx - 1)
                * profile.transition_prob(Transition::DM as usize, profile_idx - 1))
            * profile.match_prob(residue, profile_idx);
        dp_matrix.set_match(target_idx, profile_idx, m_val);

        // insert at last model position is impossible (insert_prob = 0 at M)
        let i_val = if profile_idx == profile.length {
            0.0
        } else {
            (dp_matrix.get_match(target_idx - 1, profile_idx)
                * profile.transition_prob(Transition::MI as usize, profile_idx)
                + dp_matrix.get_insert(target_idx - 1, profile_idx)
                    * profile.transition_prob(Transition::II as usize, profile_idx))
                * profile.insert_prob(residue, profile_idx)
        };
        dp_matrix.set_insert(target_idx, profile_idx, i_val);

        let d_val = dp_matrix.get_match(target_idx, profile_idx - 1)
            * profile.transition_prob(Transition::MD as usize, profile_idx - 1)
            + dp_matrix.get_delete(target_idx, profile_idx - 1)
                * profile.transition_prob(Transition::DD as usize, profile_idx - 1);
        dp_matrix.set_delete(target_idx, profile_idx, d_val);

        core_max = core_max.max(m_val).max(i_val);

        // unrolled E state
        let e_prev = dp_matrix.get_special(target_idx, Profile::E_IDX);
        dp_matrix.set_special(target_idx, Profile::E_IDX, e_prev + m_val + d_val);

        // C state
        let c_val = dp_matrix.get_special(target_idx - 1, Profile::C_IDX)
            * profile.special_transition_prob(Profile::C_IDX, Profile::SPECIAL_LOOP_IDX)
            + dp_matrix.get_special(target_idx, Profile::E_IDX)
                * profile.special_transition_prob(Profile::E_IDX, Profile::SPECIAL_MOVE_IDX);
        dp_matrix.set_special(target_idx, Profile::C_IDX, c_val);

        // N state
        let n_val = dp_matrix.get_special(target_idx - 1, Profile::N_IDX)
            * profile.special_transition_prob(Profile::N_IDX, Profile::SPECIAL_LOOP_IDX);
        dp_matrix.set_special(target_idx, Profile::N_IDX, n_val);

        // B state
        dp_matrix.set_special(
            target_idx,
            Profile::B_IDX,
            n_val * profile.special_transition_prob(Profile::N_IDX, Profile::SPECIAL_MOVE_IDX),
        );

        // Scale based on core states to keep match probability alive through
        // insertions. But if C/N are so much larger that scaling by core max
        // would push them past f32 range, fall back to the overall max.
        let special_max = dp_matrix
            .get_special(target_idx, Profile::C_IDX)
            .max(dp_matrix.get_special(target_idx, Profile::N_IDX));
        let max_val = if core_max > 0.0 && special_max / core_max < 1e30 {
            core_max
        } else {
            core_max.max(special_max)
        };

        if max_val > 0.0 {
            let inv = 1.0 / max_val;
            let lp = bounds.left_row_bounds[target_idx];
            let rp = bounds.right_row_bounds[target_idx];
            let (m_sl, i_sl, d_sl) = dp_matrix.core_slices_mut(target_idx, lp, rp);
            m_sl.iter_mut().for_each(|v| *v *= inv);
            i_sl.iter_mut().for_each(|v| *v *= inv);
            d_sl.iter_mut().for_each(|v| *v *= inv);
            dp_matrix.set_special(
                target_idx,
                Profile::N_IDX,
                dp_matrix.get_special(target_idx, Profile::N_IDX) * inv,
            );
            dp_matrix.set_special(
                target_idx,
                Profile::B_IDX,
                dp_matrix.get_special(target_idx, Profile::B_IDX) * inv,
            );
            dp_matrix.set_special(
                target_idx,
                Profile::E_IDX,
                dp_matrix.get_special(target_idx, Profile::E_IDX) * inv,
            );
            dp_matrix.set_special(
                target_idx,
                Profile::C_IDX,
                dp_matrix.get_special(target_idx, Profile::C_IDX) * inv,
            );
            xsc += max_val.ln();
        }

        // Store the per-row scale factor for use in posterior
        dp_matrix.set_row_scale(target_idx, max_val);
    }

    // Background correction for positions outside the cloud
    let aligned_target_length = bounds.seq_end - bounds.seq_start + 1;
    let unaligned_target_length = target.length - aligned_target_length;
    let background_correction = unaligned_target_length as f32
        * profile.special_transition_score(Profile::N_IDX, Profile::SPECIAL_LOOP_IDX);

    let final_c = dp_matrix.get_special(bounds.seq_end, Profile::C_IDX);
    let c_to_exit_score =
        profile.special_transition_score(Profile::C_IDX, Profile::SPECIAL_MOVE_IDX);

    Nats(final_c.ln() + xsc + background_correction + c_to_exit_score)
}
