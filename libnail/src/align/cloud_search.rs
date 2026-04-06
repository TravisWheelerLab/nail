use std::fmt::Display;

use crate::{
    align::structs::{Cloud, Seed, TraceCursor},
    alphabet::AminoAcid,
    max_f32,
    structs::{
        profile::Transition,
        Profile, Sequence,
    },
    util::MaxAssign,
};

use super::{
    structs::{Ad, BackgroundState::*, Bound, Cell, CloudMatrix, CoreState::*, NewDpMatrix},
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
pub fn compute_backward_cells<M>(
    prf: &Profile,
    seq: &Sequence,
    mx: &mut M,
    prf_idx: usize,
    seq_idx: usize,
    bg_scale: f32,
) where
    M: NewDpMatrix,
{
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

pub fn cloud_search_bwd<M>(
    prf: &Profile,
    seq: &Sequence,
    seed: &Seed,
    mx: &mut M,
    params: &CloudSearchParams,
    cloud: &mut Cloud,
) -> CloudSearchResults
where
    M: CloudMatrix,
{
    let mut num_cells = 0usize;
    let mut max_log_score = -f32::INFINITY;
    let mut max_log_score_within = -f32::INFINITY;

    let prf_end = seed.prf_end.min(prf.length);
    let seq_end = seed.seq_end.min(seq.length);

    let ad_start = prf_end + seq_end;
    let seed_ad_start = seed.seq_start + seed.prf_start;
    let min_ad = 2usize;
    let mut gamma_ad = (ad_start - params.gamma).max(min_ad);

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
    let mut trace = TraceCursor::new_backward(&seed.cigar, prf_end, seq_end);

    for idx in (min_ad..=cell_start.idx()).rev() {
        let co_located_bound = &cloud[Ad((idx + 3).min(seq.length + prf.length))];
        mx.reset_ad(co_located_bound);

        let bound = &mut cloud[Ad(idx)];
        let mut max_stored_in_ad = 0.0f32;
        bound.iter().for_each(|c| {
            compute_backward_cells(prf, seq, mx, c.prf_idx, c.seq_idx, bg_scale);
            max_stored_in_ad = max_stored_in_ad
                .max(mx[c.m_cell()])
                .max(mx[c.i_cell()])
                .max(mx[c.d_cell()]);
            num_cells += 1;
        });

        // Per-AD scaling of core cells
        if max_stored_in_ad > 0.0 {
            let inv = 1.0 / max_stored_in_ad;
            bound.iter().for_each(|c| {
                mx[c.m_cell()] *= inv;
                mx[c.i_cell()] *= inv;
                mx[c.d_cell()] *= inv;
            });
            cumulative_log_scale += max_stored_in_ad.ln();
            bg_scale = (-cumulative_log_scale).exp();
        }

        let max_log_in_ad = if max_stored_in_ad > 0.0 {
            cumulative_log_scale
        } else {
            -f32::INFINITY
        };

        let trace_on_ad = trace.active && trace.ad() == idx;
        if trace_on_ad {
            if idx <= gamma_ad {
                let prob_thresh = (max_log_in_ad - params.alpha - cumulative_log_scale).exp();
                mx.trim_ad(bound, prob_thresh);
            }

            let tp = trace.prf_idx;
            let ts = trace.seq_idx;
            if tp >= 1 && tp <= prf.length && ts >= 1 && ts <= seq.length {
                let band = 3usize;
                let left_prf = (tp + band).min(prf.length);
                let left_seq = idx.saturating_sub(left_prf).max(1);
                let left_prf = idx - left_seq;
                let right_seq = (ts + band).min(seq.length);
                let right_prf = idx.saturating_sub(right_seq).max(1);
                let right_seq = idx - right_prf;

                if bound.is_empty() {
                    bound.0 = Cell { prf_idx: left_prf, seq_idx: left_seq };
                    bound.1 = Cell { prf_idx: right_prf, seq_idx: right_seq };
                } else {
                    if left_prf > bound.0.prf_idx {
                        bound.0 = Cell { prf_idx: left_prf, seq_idx: left_seq };
                    }
                    if right_seq > bound.1.seq_idx {
                        bound.1 = Cell { prf_idx: right_prf, seq_idx: right_seq };
                    }
                }
            }

            max_log_score = max_log_in_ad;
            if idx > seed_ad_start {
                max_log_score_within.max_assign(max_log_score);
            }
            let was_active = trace.active;
            trace.advance_backward();
            if was_active && !trace.active {
                gamma_ad = idx.saturating_sub(prf.length).max(min_ad);
                max_log_score = max_log_in_ad;
            }
        } else {
            max_log_score.max_assign(max_log_in_ad);
            if idx > seed_ad_start {
                max_log_score_within.max_assign(max_log_score);
            }
            if idx <= gamma_ad {
                let thresh_log = params.trim_thresh(max_log_in_ad, max_log_score);
                let prob_thresh = (thresh_log - cumulative_log_scale).exp();
                mx.trim_ad(bound, prob_thresh);
                if bound.is_empty() {
                    break;
                }
            }
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
pub fn compute_forward_cells<M>(
    prf: &Profile,
    seq: &Sequence,
    mx: &mut M,
    prf_idx: usize,
    seq_idx: usize,
    bg_scale: f32,
) where
    M: NewDpMatrix,
{
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

pub fn cloud_search_fwd<M>(
    prf: &Profile,
    seq: &Sequence,
    seed: &Seed,
    mx: &mut M,
    params: &CloudSearchParams,
    cloud: &mut Cloud,
) -> CloudSearchResults
where
    M: CloudMatrix,
{
    let mut num_cells = 0usize;
    let mut max_log_score = -f32::INFINITY;
    let mut max_log_score_within = -f32::INFINITY;

    let ad_start = seed.seq_start + seed.prf_start;
    let seed_ad_end = seed.seq_end + seed.prf_end;
    let mut gamma_ad = ad_start + params.gamma;
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
    let mut trace = TraceCursor::new_forward(&seed.cigar, seed.prf_start, seed.seq_start);

    for idx in ad_start..=max_ad {
        let co_located_bound = &cloud[Ad(idx.saturating_sub(3))];
        mx.reset_ad(co_located_bound);

        let bound = &mut cloud[Ad(idx)];
        let mut max_stored_in_ad = 0.0f32;
        bound.iter().for_each(|c| {
            compute_forward_cells(prf, seq, mx, c.prf_idx, c.seq_idx, bg_scale);
            max_stored_in_ad = max_stored_in_ad
                .max(mx[c.m_cell()])
                .max(mx[c.i_cell()])
                .max(mx[c.d_cell()]);
            num_cells += 1;
        });

        // Per-AD scaling of core cells
        if max_stored_in_ad > 0.0 {
            let inv = 1.0 / max_stored_in_ad;
            bound.iter().for_each(|c| {
                mx[c.m_cell()] *= inv;
                mx[c.i_cell()] *= inv;
                mx[c.d_cell()] *= inv;
            });
            cumulative_log_scale += max_stored_in_ad.ln();
            bg_scale = (-cumulative_log_scale).exp();
        }

        let max_log_in_ad = if max_stored_in_ad > 0.0 {
            cumulative_log_scale
        } else {
            -f32::INFINITY
        };

        let trace_on_ad = trace.active && trace.ad() == idx;
        if trace_on_ad {
            // While the trace covers this AD: only alpha-prune (no beta),
            // enforce a minimum ±3 band around the trace, and reset max_score.
            if idx >= gamma_ad {
                let prob_thresh = (max_log_in_ad - params.alpha - cumulative_log_scale).exp();
                mx.trim_ad(bound, prob_thresh);
            }

            let tp = trace.prf_idx;
            let ts = trace.seq_idx;
            if tp >= 1 && tp <= prf.length && ts >= 1 && ts <= seq.length {
                let band = 3usize;
                let left_prf = (tp + band).min(prf.length);
                let left_seq = idx.saturating_sub(left_prf).max(1);
                let left_prf = idx - left_seq;
                let right_seq = (ts + band).min(seq.length);
                let right_prf = idx.saturating_sub(right_seq).max(1);
                let right_seq = idx - right_prf;

                if bound.is_empty() {
                    bound.0 = Cell { prf_idx: left_prf, seq_idx: left_seq };
                    bound.1 = Cell { prf_idx: right_prf, seq_idx: right_seq };
                } else {
                    if left_prf > bound.0.prf_idx {
                        bound.0 = Cell { prf_idx: left_prf, seq_idx: left_seq };
                    }
                    if right_seq > bound.1.seq_idx {
                        bound.1 = Cell { prf_idx: right_prf, seq_idx: right_seq };
                    }
                }
            }

            max_log_score = max_log_in_ad;
            if idx < seed_ad_end {
                max_log_score_within.max_assign(max_log_score);
            }
            let was_active = trace.active;
            trace.advance_forward();
            if was_active && !trace.active {
                gamma_ad = idx + prf.length;
                max_log_score = max_log_in_ad;
            }
        } else {
            max_log_score.max_assign(max_log_in_ad);
            if idx < seed_ad_end {
                max_log_score_within.max_assign(max_log_score);
            }
            if idx >= gamma_ad {
                let thresh_log = params.trim_thresh(max_log_in_ad, max_log_score);
                let prob_thresh = (thresh_log - cumulative_log_scale).exp();
                mx.trim_ad(bound, prob_thresh);
                if bound.is_empty() {
                    break;
                }
            }
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
