use crate::align::structs::{DpMatrix, RowBounds};
use crate::structs::profile::Transition;
use crate::structs::{Profile, Sequence};

/// Probability-space Backward algorithm with per-row scaling.
///
/// All DP values are stored as probabilities. Per-row scaling prevents underflow.
/// Scale factors are stored in the matrix for use by the posterior computation.
pub fn backward(
    profile: &Profile,
    target: &Sequence,
    dp_matrix: &mut impl DpMatrix,
    row_bounds: &RowBounds,
) {
    // Terminal conditions in probability space
    dp_matrix.set_special(row_bounds.seq_end, Profile::B_IDX, 0.0);
    dp_matrix.set_special(row_bounds.seq_end, Profile::N_IDX, 0.0);

    dp_matrix.set_special(
        row_bounds.seq_end,
        Profile::C_IDX,
        profile.special_transition_prob(Profile::C_IDX, Profile::SPECIAL_MOVE_IDX),
    );

    // E[T] = C->T * E->C (in prob-space: product of the two)
    dp_matrix.set_special(
        row_bounds.seq_end,
        Profile::E_IDX,
        profile.special_transition_prob(Profile::C_IDX, Profile::SPECIAL_MOVE_IDX)
            * profile.special_transition_prob(Profile::E_IDX, Profile::SPECIAL_MOVE_IDX),
    );

    let profile_start_on_last_row = row_bounds.left_row_bounds[row_bounds.seq_end];
    let profile_end_on_last_row = row_bounds.right_row_bounds[row_bounds.seq_end];

    // Last cell on last row
    dp_matrix.set_match(
        row_bounds.seq_end,
        profile_end_on_last_row,
        dp_matrix.get_special(row_bounds.seq_end, Profile::E_IDX),
    );
    dp_matrix.set_delete(
        row_bounds.seq_end,
        profile_end_on_last_row,
        dp_matrix.get_special(row_bounds.seq_end, Profile::E_IDX),
    );
    dp_matrix.set_insert(row_bounds.seq_end, profile_end_on_last_row, 0.0);

    // Fill remaining cells on last row (right to left)
    for profile_idx in (profile_start_on_last_row..profile_end_on_last_row).rev() {
        dp_matrix.set_match(
            row_bounds.seq_end,
            profile_idx,
            dp_matrix.get_special(row_bounds.seq_end, Profile::E_IDX)
                + dp_matrix.get_delete(row_bounds.seq_end, profile_idx + 1)
                    * profile.transition_prob(Transition::MD as usize, profile_idx),
        );
        dp_matrix.set_insert(row_bounds.seq_end, profile_idx, 0.0);
        dp_matrix.set_delete(
            row_bounds.seq_end,
            profile_idx,
            dp_matrix.get_special(row_bounds.seq_end, Profile::E_IDX)
                + dp_matrix.get_delete(row_bounds.seq_end, profile_idx + 1)
                    * profile.transition_prob(Transition::DD as usize, profile_idx),
        );
    }

    // Scale the last row
    {
        let mut core_max = 0.0f32;
        for p in profile_start_on_last_row..=profile_end_on_last_row {
            core_max = core_max.max(dp_matrix.get_match(row_bounds.seq_end, p));
        }
        let special_max = dp_matrix.get_special(row_bounds.seq_end, Profile::C_IDX)
            .max(dp_matrix.get_special(row_bounds.seq_end, Profile::N_IDX));
        let max_val = if core_max > 0.0 && special_max / core_max < 1e30 {
            core_max
        } else {
            core_max.max(special_max)
        };
        if max_val > 0.0 {
            let inv = 1.0 / max_val;
            for p in profile_start_on_last_row..=profile_end_on_last_row {
                dp_matrix.set_match(
                    row_bounds.seq_end,
                    p,
                    dp_matrix.get_match(row_bounds.seq_end, p) * inv,
                );
                dp_matrix.set_insert(
                    row_bounds.seq_end,
                    p,
                    dp_matrix.get_insert(row_bounds.seq_end, p) * inv,
                );
                dp_matrix.set_delete(
                    row_bounds.seq_end,
                    p,
                    dp_matrix.get_delete(row_bounds.seq_end, p) * inv,
                );
            }
            dp_matrix.set_special(
                row_bounds.seq_end,
                Profile::N_IDX,
                dp_matrix.get_special(row_bounds.seq_end, Profile::N_IDX) * inv,
            );
            dp_matrix.set_special(
                row_bounds.seq_end,
                Profile::B_IDX,
                dp_matrix.get_special(row_bounds.seq_end, Profile::B_IDX) * inv,
            );
            dp_matrix.set_special(
                row_bounds.seq_end,
                Profile::E_IDX,
                dp_matrix.get_special(row_bounds.seq_end, Profile::E_IDX) * inv,
            );
            dp_matrix.set_special(
                row_bounds.seq_end,
                Profile::C_IDX,
                dp_matrix.get_special(row_bounds.seq_end, Profile::C_IDX) * inv,
            );
        }
        dp_matrix.set_row_scale(row_bounds.seq_end, max_val);
    }

    // Main recursion (rows in reverse)
    for target_idx in (row_bounds.seq_start..row_bounds.seq_end).rev() {
        let current_residue = target.digital_bytes[target_idx + 1] as usize;
        let profile_start_on_current_row = row_bounds.left_row_bounds[target_idx];
        let profile_end_on_current_row = row_bounds.right_row_bounds[target_idx];

        // B state: sum over all profile positions of
        //   M[t+1][p] * BM_prob[p-1] * emit_M[p][residue]
        let mut b_val = dp_matrix.get_match(target_idx + 1, profile_start_on_current_row)
            * profile.transition_prob(
                Transition::BM as usize,
                profile_start_on_current_row - 1,
            )
            * profile.match_prob(current_residue, profile_start_on_current_row);

        for profile_idx in (profile_start_on_current_row + 1)..=profile_end_on_current_row {
            b_val += dp_matrix.get_match(target_idx + 1, profile_idx)
                * profile.transition_prob(Transition::BM as usize, profile_idx - 1)
                * profile.match_prob(current_residue, profile_idx);
        }
        dp_matrix.set_special(target_idx, Profile::B_IDX, b_val);

        // C state
        let c_val = dp_matrix.get_special(target_idx + 1, Profile::C_IDX)
            * profile.special_transition_prob(Profile::C_IDX, Profile::SPECIAL_LOOP_IDX);
        dp_matrix.set_special(target_idx, Profile::C_IDX, c_val);

        // E state
        let e_val = c_val
            * profile.special_transition_prob(Profile::E_IDX, Profile::SPECIAL_MOVE_IDX);
        dp_matrix.set_special(target_idx, Profile::E_IDX, e_val);

        // N state
        let n_val = dp_matrix.get_special(target_idx + 1, Profile::N_IDX)
            * profile.special_transition_prob(Profile::N_IDX, Profile::SPECIAL_LOOP_IDX)
            + b_val
                * profile.special_transition_prob(Profile::N_IDX, Profile::SPECIAL_MOVE_IDX);
        dp_matrix.set_special(target_idx, Profile::N_IDX, n_val);

        // Last profile position: M = E, I = 0, D = E
        dp_matrix.set_match(target_idx, profile_end_on_current_row, e_val);
        dp_matrix.set_insert(target_idx, profile_end_on_current_row, 0.0);
        dp_matrix.set_delete(target_idx, profile_end_on_current_row, e_val);

        // Core states (right to left)
        for profile_idx in (profile_start_on_current_row..=profile_end_on_current_row).rev() {
            let m_next_emit = dp_matrix.get_match(target_idx + 1, profile_idx + 1)
                * profile.match_prob(current_residue, profile_idx + 1);
            let i_next_emit = dp_matrix.get_insert(target_idx + 1, profile_idx)
                * profile.insert_prob(current_residue, profile_idx);

            // match state
            dp_matrix.set_match(
                target_idx,
                profile_idx,
                m_next_emit
                    * profile.transition_prob(Transition::MM as usize, profile_idx)
                    + i_next_emit
                        * profile.transition_prob(Transition::MI as usize, profile_idx)
                    + e_val
                    + dp_matrix.get_delete(target_idx, profile_idx + 1)
                        * profile.transition_prob(Transition::MD as usize, profile_idx),
            );

            // insert state
            dp_matrix.set_insert(
                target_idx,
                profile_idx,
                m_next_emit
                    * profile.transition_prob(Transition::IM as usize, profile_idx)
                    + i_next_emit
                        * profile.transition_prob(Transition::II as usize, profile_idx),
            );

            // delete state
            dp_matrix.set_delete(
                target_idx,
                profile_idx,
                m_next_emit
                    * profile.transition_prob(Transition::DM as usize, profile_idx)
                    + dp_matrix.get_delete(target_idx, profile_idx + 1)
                        * profile.transition_prob(Transition::DD as usize, profile_idx)
                    + e_val,
            );
        }

        // Scale this row — see forward.rs for rationale on core-only scaling
        let mut core_max = 0.0f32;
        for p in profile_start_on_current_row..=profile_end_on_current_row {
            core_max = core_max.max(dp_matrix.get_match(target_idx, p));
            core_max = core_max.max(dp_matrix.get_insert(target_idx, p));
        }
        let special_max = dp_matrix.get_special(target_idx, Profile::C_IDX)
            .max(dp_matrix.get_special(target_idx, Profile::N_IDX));
        let max_val = if core_max > 0.0 && special_max / core_max < 1e30 {
            core_max
        } else {
            core_max.max(special_max)
        };
        if max_val > 0.0 {
            let inv = 1.0 / max_val;
            for p in profile_start_on_current_row..=profile_end_on_current_row {
                dp_matrix.set_match(
                    target_idx,
                    p,
                    dp_matrix.get_match(target_idx, p) * inv,
                );
                dp_matrix.set_insert(
                    target_idx,
                    p,
                    dp_matrix.get_insert(target_idx, p) * inv,
                );
                dp_matrix.set_delete(
                    target_idx,
                    p,
                    dp_matrix.get_delete(target_idx, p) * inv,
                );
            }
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
        }
        dp_matrix.set_row_scale(target_idx, max_val);
    }

    // Fill row seq_start - 1 (the "before first position" row)
    let first_target_character = target.digital_bytes[row_bounds.seq_start] as usize;
    let profile_start_in_first_row = row_bounds.left_row_bounds[row_bounds.seq_start];
    let profile_end_in_first_row = row_bounds.right_row_bounds[row_bounds.seq_start];

    let mut b_val = dp_matrix.get_match(row_bounds.seq_start, profile_start_in_first_row)
        * profile.transition_prob(Transition::BM as usize, 0)
        * profile.match_prob(first_target_character, 1);

    for profile_idx in (profile_start_in_first_row + 1)..=profile_end_in_first_row {
        b_val += dp_matrix.get_match(row_bounds.seq_start, profile_idx)
            * profile.transition_prob(Transition::BM as usize, profile_idx - 1)
            * profile.match_prob(first_target_character, profile_idx);
    }
    dp_matrix.set_special(row_bounds.seq_start - 1, Profile::B_IDX, b_val);

    dp_matrix.set_special(row_bounds.seq_start - 1, Profile::C_IDX, 0.0);
    dp_matrix.set_special(row_bounds.seq_start - 1, Profile::E_IDX, 0.0);
    dp_matrix.set_special(
        row_bounds.seq_start - 1,
        Profile::N_IDX,
        dp_matrix.get_special(row_bounds.seq_start, Profile::N_IDX)
            * profile.special_transition_prob(Profile::N_IDX, Profile::SPECIAL_LOOP_IDX)
            + b_val
                * profile.special_transition_prob(Profile::N_IDX, Profile::SPECIAL_MOVE_IDX),
    );

    for profile_idx in (profile_start_in_first_row..=profile_end_in_first_row).rev() {
        dp_matrix.set_match(row_bounds.seq_start - 1, profile_idx, 0.0);
        dp_matrix.set_insert(row_bounds.seq_start - 1, profile_idx, 0.0);
        dp_matrix.set_delete(row_bounds.seq_start - 1, profile_idx, 0.0);
    }
}
