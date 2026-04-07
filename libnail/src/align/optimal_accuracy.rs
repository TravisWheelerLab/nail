use crate::align::structs::{DpMatrix, DpMatrixSparse, RowBounds};
use crate::structs::profile::Transition;
use crate::structs::Profile;

pub fn optimal_accuracy(
    profile: &Profile,
    posterior_matrix: &DpMatrixSparse,
    optimal_matrix: &mut DpMatrixSparse,
    row_bounds: &RowBounds,
) {
    // initialization of the zero row
    optimal_matrix.set_special(row_bounds.seq_start - 1, Profile::N_IDX, 0.0);
    optimal_matrix.set_special(row_bounds.seq_start - 1, Profile::B_IDX, 0.0);
    optimal_matrix.set_special(row_bounds.seq_start - 1, Profile::E_IDX, -f32::INFINITY);
    optimal_matrix.set_special(row_bounds.seq_start - 1, Profile::C_IDX, -f32::INFINITY);
    optimal_matrix.set_special(row_bounds.seq_start - 1, Profile::J_IDX, -f32::INFINITY);

    let profile_start_in_first_row = row_bounds.left_row_bounds[row_bounds.seq_start];
    let profile_end_in_first_row = row_bounds.right_row_bounds[row_bounds.seq_start];

    for profile_idx in (profile_start_in_first_row - 1)..=profile_end_in_first_row {
        optimal_matrix.set_match(row_bounds.seq_start - 1, profile_idx, -f32::INFINITY);
        optimal_matrix.set_insert(row_bounds.seq_start - 1, profile_idx, -f32::INFINITY);
        optimal_matrix.set_delete(row_bounds.seq_start - 1, profile_idx, -f32::INFINITY);
    }

    // Scratch buffers for prev-row reads (safe bounds handling).
    // Sized to cover [lp-1, rp] (length + 2 elements worst case).
    let buf_cap = profile.length + 2;
    let mut prev_m_buf = vec![0.0f32; buf_cap];
    let mut prev_i_buf = vec![0.0f32; buf_cap];
    let mut prev_d_buf = vec![0.0f32; buf_cap];

    for target_idx in row_bounds.seq_start..=row_bounds.seq_end {
        let lp = row_bounds.left_row_bounds[target_idx];
        let rp = row_bounds.right_row_bounds[target_idx];
        let len = rp - lp + 1;
        let b_prev = optimal_matrix.get_special(target_idx - 1, Profile::B_IDX);

        optimal_matrix.set_match(target_idx, lp - 1, -f32::INFINITY);
        optimal_matrix.set_insert(target_idx, lp - 1, -f32::INFINITY);
        optimal_matrix.set_delete(target_idx, lp - 1, -f32::INFINITY);

        // Fill scratch buffers: prev-row at [lp-1, rp].
        // Covers M-predecessor range [lp-1, rp-1] and I-predecessor range [lp, rp].
        let buf_len = len + 1;
        for i in 0..buf_len {
            let p = lp - 1 + i;
            prev_m_buf[i] = optimal_matrix.get_match(target_idx - 1, p);
            prev_i_buf[i] = optimal_matrix.get_insert(target_idx - 1, p);
            prev_d_buf[i] = optimal_matrix.get_delete(target_idx - 1, p);
        }

        // --- M pass ---
        // M[t][p] = max(d_MM*prev_M[p-1], d_IM*prev_I[p-1], d_DM*prev_D[p-1], d_BM*b_prev) + post_M[p]
        //
        // Factoring is valid because transition deltas are 1.0 or f32::MIN_POSITIVE, all
        // prev values and posteriors are >= 0, so max(d*(v+post)) == max(d*v) + post.
        // No within-row dependencies; auto-vectorizes.
        {
            let mm_d = profile.trans_delta_slice(Transition::MM as usize, lp - 1, rp - 1);
            let im_d = profile.trans_delta_slice(Transition::IM as usize, lp - 1, rp - 1);
            let dm_d = profile.trans_delta_slice(Transition::DM as usize, lp - 1, rp - 1);
            let bm_d = profile.trans_delta_slice(Transition::BM as usize, lp - 1, rp - 1);
            let post_m = posterior_matrix.m_slice(target_idx, lp, rp);
            let (m_sl, _, _) = optimal_matrix.core_slices_mut(target_idx, lp, rp);
            for i in 0..len {
                let max_pred = (mm_d[i] * prev_m_buf[i])
                    .max(im_d[i] * prev_i_buf[i])
                    .max(dm_d[i] * prev_d_buf[i])
                    .max(bm_d[i] * b_prev);
                m_sl[i] = max_pred + post_m[i];
            }
        }

        // --- D pass: D[t][p] = max(d_MD*M[t][p-1], d_DD*D[t][p-1])  (serial, left-to-right) ---
        // Guard cells at lp-1 were set to -inf above, so D[t][lp] = -inf.
        {
            let md_d = profile.trans_delta_slice(Transition::MD as usize, lp - 1, rp - 1);
            let dd_d = profile.trans_delta_slice(Transition::DD as usize, lp - 1, rp - 1);
            // Read guard cells before taking the mutable borrow.
            let m_guard = optimal_matrix.get_match(target_idx, lp - 1);
            let d_guard = optimal_matrix.get_delete(target_idx, lp - 1);
            let (m_sl, _, d_sl) = optimal_matrix.core_slices_mut(target_idx, lp, rp);
            d_sl[0] = (md_d[0] * m_guard).max(dd_d[0] * d_guard);
            for i in 1..len {
                d_sl[i] = (md_d[i] * m_sl[i - 1]).max(dd_d[i] * d_sl[i - 1]);
            }
        }

        // --- E = max over M[lp..=rp] and D[rp] ---
        let e_val: f32;
        {
            let (m_sl, _, d_sl) = optimal_matrix.core_slices_mut(target_idx, lp, rp);
            e_val = m_sl.iter().cloned().fold(-f32::INFINITY, f32::max)
                .max(d_sl[len - 1]);
        }
        optimal_matrix.set_special(target_idx, Profile::E_IDX, e_val);

        // --- I pass ---
        // I[t][p] = max(d_MI*prev_M[p], d_II*prev_I[p]) + post_I[p]
        // prev_{m,i}_buf[1..=len] correspond to [lp, rp] (no shift).
        // No within-row dependencies; auto-vectorizes.
        {
            let mi_d = profile.trans_delta_slice(Transition::MI as usize, lp, rp);
            let ii_d = profile.trans_delta_slice(Transition::II as usize, lp, rp);
            let post_i = posterior_matrix.i_slice(target_idx, lp, rp);
            let (_, i_sl, _) = optimal_matrix.core_slices_mut(target_idx, lp, rp);
            for i in 0..len {
                i_sl[i] = (mi_d[i] * prev_m_buf[i + 1]).max(ii_d[i] * prev_i_buf[i + 1])
                    + post_i[i];
            }
        }

        // Special states (unchanged from original)
        // a comment from hmmer:
        //   now the special states; it's important that E is already done, and B is done after N,J
        optimal_matrix.set_special(
            target_idx,
            Profile::J_IDX,
            (profile.special_transition_score_delta(Profile::J_IDX, Profile::SPECIAL_LOOP_IDX)
                * (optimal_matrix.get_special(target_idx - 1, Profile::J_IDX)
                    + posterior_matrix.get_special(target_idx, Profile::J_IDX)))
                .max(
                    profile.special_transition_score_delta(
                        Profile::E_IDX,
                        Profile::SPECIAL_LOOP_IDX,
                    ) * e_val,
                ),
        );

        optimal_matrix.set_special(
            target_idx,
            Profile::C_IDX,
            (profile.special_transition_score_delta(Profile::C_IDX, Profile::SPECIAL_LOOP_IDX)
                * (optimal_matrix.get_special(target_idx - 1, Profile::C_IDX)
                    + posterior_matrix.get_special(target_idx, Profile::C_IDX)))
                .max(
                    profile.special_transition_score_delta(
                        Profile::E_IDX,
                        Profile::SPECIAL_MOVE_IDX,
                    ) * e_val,
                ),
        );

        optimal_matrix.set_special(
            target_idx,
            Profile::N_IDX,
            profile.special_transition_score_delta(Profile::N_IDX, Profile::SPECIAL_LOOP_IDX)
                * (optimal_matrix.get_special(target_idx - 1, Profile::N_IDX)
                    + posterior_matrix.get_special(target_idx, Profile::N_IDX)),
        );

        optimal_matrix.set_special(
            target_idx,
            Profile::B_IDX,
            (profile.special_transition_score_delta(Profile::N_IDX, Profile::SPECIAL_MOVE_IDX)
                * optimal_matrix.get_special(target_idx, Profile::N_IDX))
                .max(
                    profile.special_transition_score_delta(
                        Profile::J_IDX,
                        Profile::SPECIAL_MOVE_IDX,
                    ) * optimal_matrix.get_special(target_idx, Profile::J_IDX),
                ),
        );
    }
}
