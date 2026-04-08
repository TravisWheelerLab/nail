use std::fmt::Display;

use crate::{
    align::structs::{Cloud, Seed},
    structs::{
        profile::Transition,
        Profile, Sequence,
    },
    util::MaxAssign,
};

use super::{
    structs::{Ad, AdMatrixLinear, BackgroundState::*, Bound, Cell, CoreState::*},
    Nats,
};

#[derive(Clone)]
pub struct CloudSearchParams {
    pub gamma: usize,
    pub alpha: f32,
    pub beta: f32,
}

impl Display for CloudSearchParams {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "γ{} α{} β{}", self.gamma, self.alpha, self.beta)
    }
}

impl Default for CloudSearchParams {
    fn default() -> Self {
        CloudSearchParams {
            gamma: 5,
            alpha: 12.0,
            beta: 20.0,
        }
    }
}

impl CloudSearchParams {
    pub fn trim_thresh(&self, a: f32, b: f32) -> f32 {
        (a - self.alpha).max(b - self.beta)
    }

    pub fn scale(&self, scalar: f32) -> Self {
        Self {
            gamma: (self.gamma as f32 * scalar) as usize,
            alpha: self.alpha * scalar,
            beta: self.beta * scalar,
        }
    }
}

#[derive(Default, Debug, PartialEq)]
pub struct CloudSearchResults {
    pub max_score: Nats,
    pub max_score_within: Nats,
    pub num_cells_computed: usize,
}

/// Probability-space backward cell computation for cloud search.
///
/// `bg_scale` converts precomputed background states (absolute probabilities)
/// to the same scale as the per-AD-scaled core states.
#[inline]
pub fn compute_backward_cells(
    prf: &Profile,
    seq: &Sequence,
    mx: &mut AdMatrixLinear,
    prf_idx: usize,
    seq_idx: usize,
    bg_scale: f32,
) {
    let residue = seq.digital_bytes[seq_idx + 1] as usize;

    // Cache emission probabilities (one lookup each, reused across states)
    let m_emit = prf.match_prob(residue, prf_idx + 1);
    let i_emit = prf.insert_prob(residue, prf_idx);

    // Read source cells once
    let m_src_cell = M(prf_idx + 1).cell_at(seq_idx + 1);
    let i_src_cell = I(prf_idx).cell_at(seq_idx + 1);
    let d_src_cell = D(prf_idx + 1).cell_at(seq_idx);

    let m_src_val = mx[m_src_cell];
    let i_src_val = mx[i_src_cell];
    let d_src_val = mx[d_src_cell];

    // E is a precomputed background state — scale it
    let e_scaled = mx[(E, seq_idx)] * bg_scale;

    // Pre-compute common products
    let m_src_emit = m_src_val * m_emit;
    let i_src_emit = i_src_val * i_emit;

    // Match backward
    let m_cell = M(prf_idx).cell_at(seq_idx);
    mx[m_cell] = m_src_emit * prf.transition_prob(Transition::MM as usize, prf_idx)
        + i_src_emit * prf.transition_prob(Transition::MI as usize, prf_idx)
        + d_src_val * prf.transition_prob(Transition::MD as usize, prf_idx)
        + e_scaled;

    // Insert backward
    let i_cell = I(prf_idx).cell_at(seq_idx);
    mx[i_cell] = m_src_emit * prf.transition_prob(Transition::IM as usize, prf_idx)
        + i_src_emit * prf.transition_prob(Transition::II as usize, prf_idx);

    // Delete backward
    let d_cell = D(prf_idx).cell_at(seq_idx);
    mx[d_cell] = m_src_emit * prf.transition_prob(Transition::DM as usize, prf_idx)
        + d_src_val * prf.transition_prob(Transition::DD as usize, prf_idx)
        + e_scaled;

    // B accumulation (B is accumulated in core-scaled space)
    mx[(B, seq_idx)] += m_src_val
        * prf.transition_prob(Transition::BM as usize, prf_idx)
        * m_emit;
}

pub fn cloud_search_bwd(
    prf: &Profile,
    seq: &Sequence,
    seed: &Seed,
    mx: &mut AdMatrixLinear,
    params: &CloudSearchParams,
    cloud: &mut Cloud,
) -> CloudSearchResults
{
    let mut num_cells = 0usize;
    let mut max_log_score = -f32::INFINITY;
    let mut max_log_score_within = -f32::INFINITY;

    let prf_end = seed.prf_end.min(prf.length);
    let seq_end = seed.seq_end.min(seq.length);

    let ad_start = prf_end + seq_end;
    let seed_ad_start = seed.seq_start + seed.prf_start;
    let min_ad = 2usize;
    let gamma_ad = (ad_start - params.gamma).max(min_ad);

    cloud.reuse(seq.length, prf.length);
    let cell_start = Cell {
        prf_idx: prf_end,
        seq_idx: seq_end,
    };

    cloud.append(Bound(cell_start, cell_start));

    // Precompute background (C, E) in probability space
    let c_loop_prob = prf.special_transition_prob(Profile::C_IDX, Profile::SPECIAL_LOOP_IDX);
    let c_move_prob = prf.special_transition_prob(Profile::C_IDX, Profile::SPECIAL_MOVE_IDX);
    let e_move_prob = prf.special_transition_prob(Profile::E_IDX, Profile::SPECIAL_MOVE_IDX);

    let init_seq_idx = cell_start.seq_idx + 1;
    let num_loops = seq.length as f32 - cell_start.seq_idx as f32 - 1.0;
    mx[(C, init_seq_idx)] = c_loop_prob.powf(num_loops) * c_move_prob;

    (1..=cell_start.seq_idx).rev().for_each(|seq_idx| {
        mx[(C, seq_idx)] = mx[(C, seq_idx + 1)] * c_loop_prob;
        mx[(E, seq_idx)] = mx[(C, seq_idx)] * e_move_prob;
    });

    let mut cumulative_log_scale = 0.0f32;
    let mut bg_scale = 1.0f32;
    let tb = &seed.trace_bounds;

    // Scratch buffers for slice-based inner loop.
    let buf_cap = seq.length + 2;
    let mut m_src_buf  = vec![0.0f32; buf_cap]; // M from AD idx+2 at seq_idx+1
    let mut i_src_buf  = vec![0.0f32; buf_cap]; // I from AD idx+1 at seq_idx+1
    let mut d_src_buf  = vec![0.0f32; buf_cap]; // D from AD idx+1 at seq_idx
    let mut e_buf      = vec![0.0f32; buf_cap]; // E[seq_idx] * bg_scale
    // Profile scratch
    let mut mm_buf     = vec![0.0f32; buf_cap];
    let mut mi_buf     = vec![0.0f32; buf_cap];
    let mut md_buf     = vec![0.0f32; buf_cap];
    let mut im_buf     = vec![0.0f32; buf_cap];
    let mut ii_buf     = vec![0.0f32; buf_cap];
    let mut dm_buf     = vec![0.0f32; buf_cap];
    let mut dd_buf     = vec![0.0f32; buf_cap];
    let mut bm_buf     = vec![0.0f32; buf_cap];
    let mut m_emit_buf = vec![0.0f32; buf_cap]; // match_prob at prf_idx+1, residue at seq_idx+1
    let mut i_emit_buf = vec![0.0f32; buf_cap]; // insert_prob at prf_idx, residue at seq_idx+1

    for idx in (min_ad..=cell_start.idx()).rev() {
        let co_located_idx = (idx + 3).min(seq.length + prf.length);
        {
            let co_located_bound = &cloud[Ad(co_located_idx)];
            if !co_located_bound.is_empty() {
                let (ss, se) = (co_located_bound.0.seq_idx, co_located_bound.1.seq_idx);
                if ss <= se {
                    let (m, i, d) = mx.core_slices_mut(co_located_idx, ss, se);
                    m.fill(0.0);
                    i.fill(0.0);
                    d.fill(0.0);
                }
            }
        }

        let bound = &cloud[Ad(idx)];
        if bound.is_empty() || bound.0.seq_idx > bound.1.seq_idx {
            let bound = &mut cloud[Ad(idx)];
            if let Some((prf_lo, prf_hi)) = tb.get(idx) {
                let seq_lo = idx.saturating_sub(prf_hi).max(1);
                let seq_hi = idx.saturating_sub(prf_lo).max(1);
                bound.0 = Cell { prf_idx: prf_hi, seq_idx: seq_lo };
                bound.1 = Cell { prf_idx: prf_lo, seq_idx: seq_hi };
                max_log_score = -f32::INFINITY;
            } else {
                break;
            }
            cloud.advance_reverse();
            continue;
        }

        let ss = bound.0.seq_idx;
        let se = bound.1.seq_idx;
        let ps = bound.0.prf_idx; // = idx - ss
        let len = se - ss + 1;

        // --- Phase 1: copy previous-AD source data into scratch buffers ---
        {
            let src_m = mx.m_slice(idx + 2, ss + 1, se + 1); // M at AD idx+2, seq_idx+1
            let src_i = mx.i_slice(idx + 1, ss + 1, se + 1); // I at AD idx+1, seq_idx+1
            let src_d = mx.d_slice(idx + 1, ss, se);          // D at AD idx+1, seq_idx
            let src_e = mx.bg_slice(E as usize, ss, se);
            m_src_buf[..len].copy_from_slice(src_m);
            i_src_buf[..len].copy_from_slice(src_i);
            d_src_buf[..len].copy_from_slice(src_d);
            for i in 0..len { e_buf[i] = src_e[i] * bg_scale; }
        }

        // --- Phase 2: profile gather ---
        // prf_idx = ps - i (decreasing), residue from seq_idx+1 = ss+i+1
        for i in 0..len {
            let prf_i   = ps - i;
            let residue = seq.digital_bytes[ss + i + 1] as usize;
            let t       = prf.prob_core_transitions[prf_i];
            mm_buf[i]     = t[Transition::MM as usize];
            mi_buf[i]     = t[Transition::MI as usize];
            md_buf[i]     = t[Transition::MD as usize];
            im_buf[i]     = t[Transition::IM as usize];
            ii_buf[i]     = t[Transition::II as usize];
            dm_buf[i]     = t[Transition::DM as usize];
            dd_buf[i]     = t[Transition::DD as usize];
            bm_buf[i]     = t[Transition::BM as usize];
            m_emit_buf[i] = prf.prob_emission_scores[Profile::MATCH_IDX][prf_i + 1][residue];
            i_emit_buf[i] = prf.prob_emission_scores[Profile::INSERT_IDX][prf_i][residue];
        }

        // --- Phase 3: compute M, I, D (auto-vectorizes) ---
        {
            let (m_out, i_out, d_out) = mx.core_slices_mut(idx, ss, se);
            for i in 0..len {
                let m_se = m_src_buf[i] * m_emit_buf[i];
                let i_se = i_src_buf[i] * i_emit_buf[i];
                m_out[i] = m_se * mm_buf[i] + i_se * mi_buf[i] + d_src_buf[i] * md_buf[i] + e_buf[i];
                i_out[i] = m_se * im_buf[i] + i_se * ii_buf[i];
                d_out[i] = m_se * dm_buf[i] + d_src_buf[i] * dd_buf[i] + e_buf[i];
            }
        }

        // --- Phase 4: B accumulation ---
        mx.accumulate_b(ss, &m_src_buf[..len], &bm_buf[..len], &m_emit_buf[..len]);

        // --- Phase 5: max over the AD (auto-vectorizes) ---
        let max_stored_in_ad = {
            let (m, i, d) = mx.core_slices_mut(idx, ss, se);
            m.iter().cloned().fold(0.0f32, f32::max)
                .max(i.iter().cloned().fold(0.0f32, f32::max))
                .max(d.iter().cloned().fold(0.0f32, f32::max))
        };
        num_cells += len;

        // Per-AD scaling of core cells
        if max_stored_in_ad > 0.0 {
            let inv = 1.0 / max_stored_in_ad;
            let (m, i, d) = mx.core_slices_mut(idx, ss, se);
            m.iter_mut().for_each(|v| *v *= inv);
            i.iter_mut().for_each(|v| *v *= inv);
            d.iter_mut().for_each(|v| *v *= inv);
            cumulative_log_scale += max_stored_in_ad.ln();
            bg_scale = (-cumulative_log_scale).exp();
        }

        let max_log_in_ad = if max_stored_in_ad > 0.0 {
            cumulative_log_scale
        } else {
            -f32::INFINITY
        };

        let bound = &mut cloud[Ad(idx)];
        if idx <= gamma_ad {
            let thresh_log = params.trim_thresh(max_log_in_ad, max_log_score);
            let prob_thresh = (thresh_log - cumulative_log_scale).exp();
            mx.trim_ad_linear(bound, prob_thresh);
        }

        if let Some((prf_lo, prf_hi)) = tb.get(idx) {
            let seq_lo = idx.saturating_sub(prf_hi).max(1);
            let seq_hi = idx.saturating_sub(prf_lo).max(1);
            if bound.is_empty() {
                bound.0 = Cell { prf_idx: prf_hi, seq_idx: seq_lo };
                bound.1 = Cell { prf_idx: prf_lo, seq_idx: seq_hi };
            } else {
                if prf_hi > bound.0.prf_idx {
                    bound.0 = Cell { prf_idx: prf_hi, seq_idx: seq_lo };
                }
                if seq_hi > bound.1.seq_idx {
                    bound.1 = Cell { prf_idx: prf_lo, seq_idx: seq_hi };
                }
            }
            max_log_score = max_log_in_ad;
        } else if bound.is_empty() {
            break;
        }

        max_log_score.max_assign(max_log_in_ad);
        if idx > seed_ad_start {
            max_log_score_within.max_assign(max_log_score);
        }

        cloud.advance_reverse();
    }

    cloud.ad_start += 1;

    // Post-cloud N state computation in probability space
    let n_loop_prob = prf.special_transition_prob(Profile::N_IDX, Profile::SPECIAL_LOOP_IDX);
    let n_move_prob = prf.special_transition_prob(Profile::N_IDX, Profile::SPECIAL_MOVE_IDX);
    (1..=cell_start.seq_idx).rev().for_each(|seq_idx| {
        mx[(N, seq_idx)] =
            mx[(N, seq_idx + 1)] * n_loop_prob + mx[(B, seq_idx)] * n_move_prob;
    });

    CloudSearchResults {
        max_score: Nats(max_log_score),
        max_score_within: Nats(max_log_score_within),
        num_cells_computed: num_cells,
    }
}

/// Probability-space forward cell computation for cloud search.
///
/// `bg_scale` converts precomputed background states (absolute probabilities)
/// to the same scale as the per-AD-scaled core states.
#[inline]
pub fn compute_forward_cells(
    prf: &Profile,
    seq: &Sequence,
    mx: &mut AdMatrixLinear,
    prf_idx: usize,
    seq_idx: usize,
    bg_scale: f32,
) {
    let residue = seq.digital_bytes[seq_idx] as usize;

    // Cache source cell values
    let m_prev = mx[(M(prf_idx - 1), seq_idx - 1)];
    let i_prev = mx[(I(prf_idx - 1), seq_idx - 1)];
    let d_prev = mx[(D(prf_idx - 1), seq_idx - 1)];

    // B is a precomputed background state — scale it
    let b_scaled = mx[(B, seq_idx - 1)] * bg_scale;

    // Match state
    let m_cell = (M(prf_idx), seq_idx);
    let m_emit = prf.match_prob(residue, prf_idx);
    mx[m_cell] = (m_prev * prf.transition_prob(Transition::MM as usize, prf_idx - 1)
        + i_prev * prf.transition_prob(Transition::IM as usize, prf_idx - 1)
        + d_prev * prf.transition_prob(Transition::DM as usize, prf_idx - 1)
        + b_scaled * prf.transition_prob(Transition::BM as usize, prf_idx - 1))
        * m_emit;

    // Insert state
    let i_cell = (I(prf_idx), seq_idx);
    let m_same = mx[(M(prf_idx), seq_idx - 1)];
    let i_same = mx[(I(prf_idx), seq_idx - 1)];
    mx[i_cell] = (m_same * prf.transition_prob(Transition::MI as usize, prf_idx)
        + i_same * prf.transition_prob(Transition::II as usize, prf_idx))
        * prf.insert_prob(residue, prf_idx);

    // Delete state
    let d_cell = (D(prf_idx), seq_idx);
    let m_left = mx[(M(prf_idx - 1), seq_idx)];
    let d_left = mx[(D(prf_idx - 1), seq_idx)];
    mx[d_cell] = m_left * prf.transition_prob(Transition::MD as usize, prf_idx - 1)
        + d_left * prf.transition_prob(Transition::DD as usize, prf_idx - 1);

    // E state accumulation (stays in core-scaled space)
    mx[(E, seq_idx)] += mx[m_cell] + mx[d_cell];
}

pub fn cloud_search_fwd(
    prf: &Profile,
    seq: &Sequence,
    seed: &Seed,
    mx: &mut AdMatrixLinear,
    params: &CloudSearchParams,
    cloud: &mut Cloud,
) -> CloudSearchResults
{
    let mut num_cells = 0usize;
    let mut max_log_score = -f32::INFINITY;
    let mut max_log_score_within = -f32::INFINITY;

    let ad_start = seed.seq_start + seed.prf_start;
    let seed_ad_end = seed.seq_end + seed.prf_end;
    let gamma_ad = ad_start + params.gamma;
    let max_ad = seq.length + prf.length;

    let cell_start = Cell {
        prf_idx: seed.prf_start,
        seq_idx: seed.seq_start,
    };

    cloud.append(Bound(cell_start, cell_start));

    // Precompute background (N, B) in probability space
    let n_loop_prob = prf.special_transition_prob(Profile::N_IDX, Profile::SPECIAL_LOOP_IDX);
    let n_move_prob = prf.special_transition_prob(Profile::N_IDX, Profile::SPECIAL_MOVE_IDX);

    let init_seq_idx = cell_start.seq_idx - 1;
    let num_loops = init_seq_idx;
    mx[(N, init_seq_idx)] = n_loop_prob.powi(num_loops as i32);
    mx[(B, init_seq_idx)] = mx[(N, init_seq_idx)] * n_move_prob;
    (cell_start.seq_idx..=seq.length).for_each(|seq_idx| {
        mx[(N, seq_idx)] = mx[(N, seq_idx - 1)] * n_loop_prob;
        mx[(B, seq_idx)] = mx[(N, seq_idx)] * n_move_prob;
    });

    let mut cumulative_log_scale = 0.0f32;
    let mut bg_scale = 1.0f32;
    let tb = &seed.trace_bounds;

    // Scratch buffers: previous-AD source data + profile gather.
    // Sized to cover the widest possible AD (seq_len + 1 elements).
    let buf_cap = seq.length + 2;
    let mut m2_buf   = vec![0.0f32; buf_cap]; // M from AD idx-2 at seq_idx-1
    let mut i2_buf   = vec![0.0f32; buf_cap]; // I from AD idx-2 at seq_idx-1
    let mut d2_buf   = vec![0.0f32; buf_cap]; // D from AD idx-2 at seq_idx-1
    let mut m1s_buf  = vec![0.0f32; buf_cap]; // M from AD idx-1 at seq_idx-1 (for I pass)
    let mut i1s_buf  = vec![0.0f32; buf_cap]; // I from AD idx-1 at seq_idx-1
    let mut m1_buf   = vec![0.0f32; buf_cap]; // M from AD idx-1 at seq_idx   (for D pass)
    let mut d1_buf   = vec![0.0f32; buf_cap]; // D from AD idx-1 at seq_idx
    let mut b_buf    = vec![0.0f32; buf_cap]; // B[seq_idx-1] * bg_scale
    // Profile scratch (indexed by position-in-AD = seq_idx - ss)
    let mut mm_buf     = vec![0.0f32; buf_cap];
    let mut im_buf     = vec![0.0f32; buf_cap];
    let mut dm_buf     = vec![0.0f32; buf_cap];
    let mut bm_buf     = vec![0.0f32; buf_cap];
    let mut mi_buf     = vec![0.0f32; buf_cap];
    let mut ii_buf     = vec![0.0f32; buf_cap];
    let mut md_buf     = vec![0.0f32; buf_cap];
    let mut dd_buf     = vec![0.0f32; buf_cap];
    let mut m_emit_buf = vec![0.0f32; buf_cap];
    let mut i_emit_buf = vec![0.0f32; buf_cap];

    for idx in ad_start..=max_ad {
        let co_located_idx = idx.saturating_sub(3);
        {
            let co_located_bound = &cloud[Ad(co_located_idx)];
            if !co_located_bound.is_empty() {
                let (ss, se) = (co_located_bound.0.seq_idx, co_located_bound.1.seq_idx);
                if ss <= se {
                    let (m, i, d) = mx.core_slices_mut(co_located_idx, ss, se);
                    m.fill(0.0);
                    i.fill(0.0);
                    d.fill(0.0);
                }
            }
        }

        let bound = &cloud[Ad(idx)];
        if bound.is_empty() || bound.0.seq_idx > bound.1.seq_idx {
            // No cells — still need trace-bound expansion and break logic below.
            let bound = &mut cloud[Ad(idx)];
            if let Some((prf_lo, prf_hi)) = tb.get(idx) {
                let seq_lo = idx.saturating_sub(prf_hi).max(1);
                let seq_hi = idx.saturating_sub(prf_lo).max(1);
                bound.0 = Cell { prf_idx: prf_hi, seq_idx: seq_lo };
                bound.1 = Cell { prf_idx: prf_lo, seq_idx: seq_hi };
                max_log_score = -f32::INFINITY;
            } else {
                break;
            }
            if cloud.ad_end + 1 >= cloud.bounds.len() { break; }
            cloud.advance_forward();
            continue;
        }

        let ss = bound.0.seq_idx; // low seq_idx (high prf_idx)
        let se = bound.1.seq_idx; // high seq_idx (low prf_idx)
        let ps = bound.0.prf_idx; // = idx - ss
        let len = se - ss + 1;

        // --- Phase 1: copy previous-AD data into scratch buffers ---
        // These are immutable borrows; we release them before the mutable write phase.
        {
            let src_m2  = mx.m_slice(idx - 2, ss - 1, se - 1);
            let src_i2  = mx.i_slice(idx - 2, ss - 1, se - 1);
            let src_d2  = mx.d_slice(idx - 2, ss - 1, se - 1);
            let src_m1s = mx.m_slice(idx - 1, ss - 1, se - 1);
            let src_i1s = mx.i_slice(idx - 1, ss - 1, se - 1);
            let src_m1  = mx.m_slice(idx - 1, ss, se);
            let src_d1  = mx.d_slice(idx - 1, ss, se);
            let src_b   = mx.bg_slice(B as usize, ss - 1, se - 1);
            m2_buf [..len].copy_from_slice(src_m2);
            i2_buf [..len].copy_from_slice(src_i2);
            d2_buf [..len].copy_from_slice(src_d2);
            m1s_buf[..len].copy_from_slice(src_m1s);
            i1s_buf[..len].copy_from_slice(src_i1s);
            m1_buf [..len].copy_from_slice(src_m1);
            d1_buf [..len].copy_from_slice(src_d1);
            // Pre-scale B
            for i in 0..len { b_buf[i] = src_b[i] * bg_scale; }
        }

        // --- Phase 2: profile gather (scatter → sequential) ---
        // For position i in AD: prf_idx = ps - i (decreasing), seq_idx = ss + i.
        for i in 0..len {
            let prf_i   = ps - i;
            let residue = seq.digital_bytes[ss + i] as usize;
            let t       = prf.prob_core_transitions[prf_i - 1]; // copy [f32; 8]
            mm_buf[i]     = t[Transition::MM as usize];
            im_buf[i]     = t[Transition::IM as usize];
            dm_buf[i]     = t[Transition::DM as usize];
            bm_buf[i]     = t[Transition::BM as usize];
            md_buf[i]     = t[Transition::MD as usize];
            dd_buf[i]     = t[Transition::DD as usize];
            let ti        = prf.prob_core_transitions[prf_i];
            mi_buf[i]     = ti[Transition::MI as usize];
            ii_buf[i]     = ti[Transition::II as usize];
            m_emit_buf[i] = prf.prob_emission_scores[Profile::MATCH_IDX][prf_i][residue];
            i_emit_buf[i] = prf.prob_emission_scores[Profile::INSERT_IDX][prf_i][residue];
        }

        // --- Phase 3: compute M, I, D (auto-vectorizes) ---
        {
            let (m_out, i_out, d_out) = mx.core_slices_mut(idx, ss, se);
            for i in 0..len {
                m_out[i] = (m2_buf[i]*mm_buf[i]
                    + i2_buf[i]*im_buf[i]
                    + d2_buf[i]*dm_buf[i]
                    + b_buf[i]*bm_buf[i])
                    * m_emit_buf[i];
            }
            for i in 0..len {
                i_out[i] = (m1s_buf[i]*mi_buf[i] + i1s_buf[i]*ii_buf[i]) * i_emit_buf[i];
            }
            for i in 0..len {
                d_out[i] = m1_buf[i]*md_buf[i] + d1_buf[i]*dd_buf[i];
            }
        }

        // --- Phase 4: E accumulation (internal method accesses both fields) ---
        mx.accumulate_e(idx, ss, se);

        // --- Phase 5: max over the AD (auto-vectorizes) ---
        let max_stored_in_ad = {
            let (m, i, d) = mx.core_slices_mut(idx, ss, se);
            m.iter().cloned().fold(0.0f32, f32::max)
                .max(i.iter().cloned().fold(0.0f32, f32::max))
                .max(d.iter().cloned().fold(0.0f32, f32::max))
        };
        num_cells += len;

        // Per-AD scaling of core cells
        if max_stored_in_ad > 0.0 {
            let inv = 1.0 / max_stored_in_ad;
            let (m, i, d) = mx.core_slices_mut(idx, ss, se);
            m.iter_mut().for_each(|v| *v *= inv);
            i.iter_mut().for_each(|v| *v *= inv);
            d.iter_mut().for_each(|v| *v *= inv);
            cumulative_log_scale += max_stored_in_ad.ln();
            bg_scale = (-cumulative_log_scale).exp();
        }

        let max_log_in_ad = if max_stored_in_ad > 0.0 {
            cumulative_log_scale
        } else {
            -f32::INFINITY
        };

        let bound = &mut cloud[Ad(idx)];
        if idx >= gamma_ad {
            let thresh_log = params.trim_thresh(max_log_in_ad, max_log_score);
            let prob_thresh = (thresh_log - cumulative_log_scale).exp();
            mx.trim_ad_linear(bound, prob_thresh);
        }

        // If the trace covers this AD, ensure the bound includes the
        // trace range and reset max_score so beta tracks current level.
        if let Some((prf_lo, prf_hi)) = tb.get(idx) {
            let seq_lo = idx.saturating_sub(prf_hi).max(1);
            let seq_hi = idx.saturating_sub(prf_lo).max(1);
            if bound.is_empty() {
                bound.0 = Cell { prf_idx: prf_hi, seq_idx: seq_lo };
                bound.1 = Cell { prf_idx: prf_lo, seq_idx: seq_hi };
            } else {
                if prf_hi > bound.0.prf_idx {
                    bound.0 = Cell { prf_idx: prf_hi, seq_idx: seq_lo };
                }
                if seq_hi > bound.1.seq_idx {
                    bound.1 = Cell { prf_idx: prf_lo, seq_idx: seq_hi };
                }
            }
            max_log_score = max_log_in_ad;
        } else if bound.is_empty() {
            break;
        }

        max_log_score.max_assign(max_log_in_ad);
        if idx < seed_ad_end {
            max_log_score_within.max_assign(max_log_score);
        }

        if cloud.ad_end + 1 >= cloud.bounds.len() {
            break;
        }
        cloud.advance_forward();
    }

    cloud.ad_end -= 1;

    // Post-cloud C state computation in probability space
    let c_loop_prob = prf.special_transition_prob(Profile::N_IDX, Profile::SPECIAL_LOOP_IDX);
    let e_move_prob = prf.special_transition_prob(Profile::E_IDX, Profile::SPECIAL_MOVE_IDX);
    (cell_start.seq_idx..=seq.length).for_each(|seq_idx| {
        mx[(C, seq_idx)] =
            mx[(C, seq_idx - 1)] * c_loop_prob + mx[(E, seq_idx)] * e_move_prob;
    });

    CloudSearchResults {
        max_score: Nats(max_log_score),
        max_score_within: Nats(max_log_score_within),
        num_cells_computed: num_cells,
    }
}
