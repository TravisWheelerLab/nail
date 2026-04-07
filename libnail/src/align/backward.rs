use crate::align::structs::{DpMatrix, DpMatrixSparse, RowBounds};
use crate::structs::profile::Transition;
use crate::structs::{Profile, Sequence};

/// Probability-space Backward algorithm with per-row scaling.
///
/// All DP values are stored as probabilities. Per-row scaling prevents underflow.
/// Scale factors are stored in the matrix for use by the posterior computation.
pub fn backward(
    profile: &Profile,
    target: &Sequence,
    dp_matrix: &mut DpMatrixSparse,
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

    // Scale the last row: I = 0 on this row, so scan only M values for core_max.
    {
        let core_max = dp_matrix
            .m_slice(row_bounds.seq_end, profile_start_on_last_row, profile_end_on_last_row)
            .iter()
            .cloned()
            .fold(0.0f32, f32::max);
        let special_max = dp_matrix.get_special(row_bounds.seq_end, Profile::C_IDX)
            .max(dp_matrix.get_special(row_bounds.seq_end, Profile::N_IDX));
        let max_val = if core_max > 0.0 && special_max / core_max < 1e30 {
            core_max
        } else {
            core_max.max(special_max)
        };
        if max_val > 0.0 {
            let inv = 1.0 / max_val;
            let (m_sl, i_sl, d_sl) = dp_matrix.core_slices_mut(
                row_bounds.seq_end,
                profile_start_on_last_row,
                profile_end_on_last_row,
            );
            m_sl.iter_mut().for_each(|v| *v *= inv);
            i_sl.iter_mut().for_each(|v| *v *= inv);
            d_sl.iter_mut().for_each(|v| *v *= inv);
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

    // Scratch buffers for next-row values, precomputed per row.
    // m_next_emit_buf[i] = M[t+1][lp+i+1] * match_emit[residue][lp+i+1]
    // i_next_emit_buf[i] = I[t+1][lp+i]   * insert_emit[residue][lp+i]
    // Sized to profile.length + 2 to cover the largest possible row.
    let buf_cap = profile.length + 2;
    let mut m_next_emit_buf = vec![0.0f32; buf_cap];
    let mut i_next_emit_buf = vec![0.0f32; buf_cap];

    // Main recursion (rows in reverse)
    for target_idx in (row_bounds.seq_start..row_bounds.seq_end).rev() {
        let current_residue = target.digital_bytes[target_idx + 1] as usize;
        let lp = row_bounds.left_row_bounds[target_idx];
        let rp = row_bounds.right_row_bounds[target_idx];
        let len = rp - lp + 1;

        // Fill scratch buffers (scalar — next-row bounds may differ so can't use slices directly).
        for i in 0..len {
            let p = lp + i;
            m_next_emit_buf[i] = dp_matrix.get_match(target_idx + 1, p + 1)
                * profile.match_prob(current_residue, p + 1);
            i_next_emit_buf[i] = dp_matrix.get_insert(target_idx + 1, p)
                * profile.insert_prob(current_residue, p);
        }

        // B state: sum over p of M[t+1][p] * BM[p-1] * match_emit[residue][p]
        // = sum of m_next_emit_buf[i] * BM[lp+i] / match_emit[lp+i+1] * match_emit[lp+i]
        // Easier to recompute directly since scratch buf uses p+1 offsets.
        let mut b_val = 0.0f32;
        for i in 0..len {
            let p = lp + i;
            b_val += dp_matrix.get_match(target_idx + 1, p)
                * profile.transition_prob(Transition::BM as usize, p - 1)
                * profile.match_prob(current_residue, p);
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

        // --- D pass: serial right-to-left (D[t][p] depends on D[t][p+1])
        // D[t][p] = m_next_emit_buf[p-lp] * DM[p] + D[t][p+1] * DD[p] + e_val
        // Last element: D[t][rp+1] = 0 (out-of-bounds guard cell)
        {
            let dm = profile.trans_slice(Transition::DM as usize, lp, rp);
            let dd = profile.trans_slice(Transition::DD as usize, lp, rp);
            let (_, _, d_sl) = dp_matrix.core_slices_mut(target_idx, lp, rp);
            let last = len - 1;
            d_sl[last] = m_next_emit_buf[last] * dm[last] + e_val; // D[t][rp+1] = 0
            for i in (0..last).rev() {
                d_sl[i] = m_next_emit_buf[i] * dm[i] + d_sl[i + 1] * dd[i] + e_val;
            }
        }

        // --- M pass: after D is computed (reads D[t][p+1] from current row)
        // M[t][p] = m_next_emit_buf[p-lp]*MM[p] + i_next_emit_buf[p-lp]*MI[p] + e_val + D[t][p+1]*MD[p]
        // Last element: D[t][rp+1] = 0, so no MD term.
        // Loop for 0..last is independent across iterations — auto-vectorizes.
        {
            let mm = profile.trans_slice(Transition::MM as usize, lp, rp);
            let mi = profile.trans_slice(Transition::MI as usize, lp, rp);
            let md = profile.trans_slice(Transition::MD as usize, lp, rp);
            let (m_sl, _, d_sl) = dp_matrix.core_slices_mut(target_idx, lp, rp);
            let last = len - 1;
            m_sl[last] = m_next_emit_buf[last] * mm[last]
                + i_next_emit_buf[last] * mi[last]
                + e_val;
            for i in 0..last {
                m_sl[i] = m_next_emit_buf[i] * mm[i]
                    + i_next_emit_buf[i] * mi[i]
                    + e_val
                    + d_sl[i + 1] * md[i];
            }
        }

        // --- I pass: fully independent across positions — auto-vectorizes.
        // I[t][p] = m_next_emit_buf[p-lp]*IM[p] + i_next_emit_buf[p-lp]*II[p]
        {
            let im = profile.trans_slice(Transition::IM as usize, lp, rp);
            let ii = profile.trans_slice(Transition::II as usize, lp, rp);
            let (_, i_sl, _) = dp_matrix.core_slices_mut(target_idx, lp, rp);
            for i in 0..len {
                i_sl[i] = m_next_emit_buf[i] * im[i] + i_next_emit_buf[i] * ii[i];
            }
        }

        // core_max for scaling
        let core_max: f32;
        {
            let (m_sl, i_sl, _) = dp_matrix.core_slices_mut(target_idx, lp, rp);
            core_max = m_sl.iter().cloned().fold(0.0f32, f32::max)
                .max(i_sl.iter().cloned().fold(0.0f32, f32::max));
        }

        // Scale this row — see forward.rs for rationale on core-only scaling
        let special_max = dp_matrix.get_special(target_idx, Profile::C_IDX)
            .max(dp_matrix.get_special(target_idx, Profile::N_IDX));
        let max_val = if core_max > 0.0 && special_max / core_max < 1e30 {
            core_max
        } else {
            core_max.max(special_max)
        };
        if max_val > 0.0 {
            let inv = 1.0 / max_val;
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
