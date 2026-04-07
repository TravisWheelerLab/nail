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

    // Scratch buffers for prev-row reads (safe bounds handling via get_match/get_insert/get_delete).
    // Sized to cover [lp-1, rp] in the worst case (length + 2 elements).
    let buf_cap = profile.length + 2;
    let mut prev_m_buf = vec![0.0f32; buf_cap];
    let mut prev_i_buf = vec![0.0f32; buf_cap];
    let mut prev_d_buf = vec![0.0f32; buf_cap];

    for target_idx in bounds.seq_start..=bounds.seq_end {
        let residue = target.digital_bytes[target_idx] as usize;
        let lp = bounds.left_row_bounds[target_idx];
        let rp = bounds.right_row_bounds[target_idx];
        let len = rp - lp + 1;
        let b_prev = dp_matrix.get_special(target_idx - 1, Profile::B_IDX);

        dp_matrix.set_special(target_idx, Profile::E_IDX, 0.0);

        // Fill scratch buffers with prev-row values at positions [lp-1, rp].
        // Covers both the M-predecessor range [lp-1, rp-1] and the I-predecessor range [lp, rp].
        // get_match/get_insert/get_delete return 0.0 for out-of-bounds accesses.
        let buf_len = len + 1; // covers [lp-1, rp]
        for i in 0..buf_len {
            let p = lp - 1 + i;
            prev_m_buf[i] = dp_matrix.get_match(target_idx - 1, p);
            prev_i_buf[i] = dp_matrix.get_insert(target_idx - 1, p);
            prev_d_buf[i] = dp_matrix.get_delete(target_idx - 1, p);
        }

        // --- M pass: M[t][p] = (prev_M[p-1]*MM + prev_I[p-1]*IM + prev_D[p-1]*DM + b_prev*BM) * emit_M[p]
        // prev_{m,i,d}_buf[0..len] correspond to positions [lp-1, rp-1] (shifted by -1 from [lp,rp]).
        // No intra-row dependencies; auto-vectorizes.
        {
            let mm = profile.trans_slice(Transition::MM as usize, lp - 1, rp - 1);
            let im = profile.trans_slice(Transition::IM as usize, lp - 1, rp - 1);
            let dm = profile.trans_slice(Transition::DM as usize, lp - 1, rp - 1);
            let bm = profile.trans_slice(Transition::BM as usize, lp - 1, rp - 1);
            let me = profile.match_emit_slice(residue, lp, rp);
            let (m_sl, _, _) = dp_matrix.core_slices_mut(target_idx, lp, rp);
            for i in 0..len {
                m_sl[i] = (prev_m_buf[i] * mm[i]
                    + prev_i_buf[i] * im[i]
                    + prev_d_buf[i] * dm[i]
                    + b_prev * bm[i])
                    * me[i];
            }
        }

        // --- I pass: I[t][p] = (prev_M[p]*MI + prev_I[p]*II) * emit_I[p]
        // prev_{m,i}_buf[1..=len] correspond to positions [lp, rp] (no shift).
        // Insert emission at profile.length is 0.0, so no special-case needed.
        // No intra-row dependencies; auto-vectorizes.
        {
            let mi = profile.trans_slice(Transition::MI as usize, lp, rp);
            let ii = profile.trans_slice(Transition::II as usize, lp, rp);
            let ie = profile.insert_emit_slice(residue, lp, rp);
            let (_, i_sl, _) = dp_matrix.core_slices_mut(target_idx, lp, rp);
            for i in 0..len {
                i_sl[i] = (prev_m_buf[i + 1] * mi[i] + prev_i_buf[i + 1] * ii[i]) * ie[i];
            }
        }

        // --- D pass: D[t][p] = M[t][p-1]*MD[p-1] + D[t][p-1]*DD[p-1]  (serial, left-to-right)
        // D[t][lp] = 0 since both guard cells (M[t][lp-1] and D[t][lp-1]) are 0.
        {
            let md = profile.trans_slice(Transition::MD as usize, lp - 1, rp - 1);
            let dd = profile.trans_slice(Transition::DD as usize, lp - 1, rp - 1);
            let (m_sl, _, d_sl) = dp_matrix.core_slices_mut(target_idx, lp, rp);
            d_sl[0] = 0.0;
            for i in 1..len {
                d_sl[i] = m_sl[i - 1] * md[i] + d_sl[i - 1] * dd[i];
            }
        }

        // --- E = sum(M + D), core_max = max(M, I) ---
        let e_val: f32;
        let core_max: f32;
        {
            let (m_sl, i_sl, d_sl) = dp_matrix.core_slices_mut(target_idx, lp, rp);
            e_val = m_sl.iter().sum::<f32>() + d_sl.iter().sum::<f32>();
            core_max = m_sl.iter().cloned().fold(0.0f32, f32::max)
                .max(i_sl.iter().cloned().fold(0.0f32, f32::max));
        }
        dp_matrix.set_special(target_idx, Profile::E_IDX, e_val);

        // C state
        let c_val = dp_matrix.get_special(target_idx - 1, Profile::C_IDX)
            * profile.special_transition_prob(Profile::C_IDX, Profile::SPECIAL_LOOP_IDX)
            + e_val * profile.special_transition_prob(Profile::E_IDX, Profile::SPECIAL_MOVE_IDX);
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
