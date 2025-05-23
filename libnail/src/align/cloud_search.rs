use crate::align::structs::{Cloud, Seed};
use crate::log_sum;
use crate::max_f32;
use crate::structs::profile::{AminoAcid, BackgroundLoop, CoreToCore, Emission};
use crate::structs::{Profile, Sequence};
use crate::util::{log_add, MaxAssign};

use super::structs::{
    Ad, BackgroundState::*, Bound, Cell, CoreState::*, NewCloudMatrix, NewDpMatrix,
};
use super::Nats;

#[derive(Clone)]
pub struct CloudSearchParams {
    pub gamma: usize,
    pub alpha: f32,
    pub beta: f32,
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
}

#[derive(Default, Debug, PartialEq)]
pub struct CloudSearchResults {
    pub max_score: Nats,
    pub max_score_within: Nats,
    pub num_cells_computed: usize,
}

#[inline]
pub fn compute_backward_cells<M>(
    prf: &Profile,
    seq: &Sequence,
    mx: &mut M,
    prf_idx: usize,
    seq_idx: usize,
) where
    M: NewDpMatrix,
{
    let residue = AminoAcid::from_u8(seq.digital_bytes[seq_idx + 1]);

    // note: the Backward recurrence is somewhat counterintuitive,
    //       since it is derived by inverting all of the transitions
    //       in the profile HMM
    //
    //       consequently, when any specific state contributes score
    //       to a Backward cell, the source cell is always the same
    //
    //       i.e., the following source states/cells are used for each
    //       of the following Backward cell computations

    let m_src = M(prf_idx + 1);
    let i_src = I(prf_idx);
    let d_src = D(prf_idx + 1);

    let m_src_cell = m_src.cell_at(seq_idx + 1);
    let i_src_cell = i_src.cell_at(seq_idx + 1);
    let d_src_cell = d_src.cell_at(seq_idx);

    let m_dest = M(prf_idx);
    let m_cell = m_dest.cell_at(seq_idx);

    mx[m_cell] = log_sum!(
        mx[m_src_cell] + prf[CoreToCore(m_dest, m_src)] + prf[Emission(m_src, residue)],
        mx[i_src_cell] + prf[CoreToCore(m_dest, i_src)] + prf[Emission(i_src, residue)],
        mx[d_src_cell] + prf[CoreToCore(m_dest, d_src)],
        mx[(E, seq_idx)]
    );

    let i_dest = I(prf_idx);
    let i_cell = i_dest.cell_at(seq_idx);

    mx[i_cell] = log_sum!(
        mx[m_src_cell] + prf[CoreToCore(i_dest, m_src)] + prf[Emission(m_src, residue)],
        mx[i_src_cell] + prf[CoreToCore(i_dest, i_src)] + prf[Emission(i_src, residue)]
    );

    let d_dest = D(prf_idx);
    let d_cell = d_dest.cell_at(seq_idx);

    mx[d_cell] = log_sum!(
        mx[m_src_cell] + prf[CoreToCore(d_dest, m_src)] + prf[Emission(m_src, residue)],
        mx[d_src_cell] + prf[CoreToCore(d_dest, d_src)],
        mx[(E, seq_idx)]
    );

    mx[(B, seq_idx)] = log_sum!(
        mx[(B, seq_idx)],
        mx[m_src_cell]
            + prf.transition_score(Profile::B_M_IDX, prf_idx)
            + prf[Emission(m_src, residue)]
    );
}

pub fn cloud_search_backward<M>(
    prf: &Profile,
    seq: &Sequence,
    seed: &Seed,
    mx: &mut M,
    params: &CloudSearchParams,
    cloud: &mut Cloud,
) -> CloudSearchResults
where
    M: NewCloudMatrix,
{
    let mut num_cells = 0usize;
    let mut max_score = -f32::INFINITY;
    let mut max_score_within = -f32::INFINITY;

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

    // precompute background (C) scores
    let seq_idx = cell_start.seq_idx + 1;
    let num_loops = seq.length as f32 - cell_start.seq_idx as f32 - 1.0;
    mx[(C, seq_idx)] = num_loops * prf[BackgroundLoop(C)]
        + prf.special_transitions[Profile::C_IDX][Profile::SPECIAL_MOVE_IDX];

    (1..=cell_start.seq_idx).rev().for_each(|seq_idx| {
        mx[(C, seq_idx)] = mx[(C, seq_idx + 1)] + prf[BackgroundLoop(C)];
        mx[(E, seq_idx)] =
            mx[(C, seq_idx)] + prf.special_transitions[Profile::E_IDX][Profile::SPECIAL_MOVE_IDX];
    });

    for idx in (min_ad..=cell_start.idx()).rev() {
        let co_located_bound = &cloud[Ad((idx + 3).min(seq.length + prf.length))];
        mx.reset_ad(co_located_bound);

        let bound = &mut cloud[Ad(idx)];
        let mut max_score_in_ad = -f32::INFINITY;
        bound.iter().for_each(|c| {
            compute_backward_cells(prf, seq, mx, c.prf_idx, c.seq_idx);
            max_score_in_ad = max_f32!(
                max_score_in_ad,
                mx[c.m_cell()],
                mx[c.i_cell()],
                mx[c.d_cell()]
            );
            num_cells += 1;
        });

        max_score.max_assign(max_score_in_ad);
        if idx > seed_ad_start {
            max_score_within.max_assign(max_score)
        }

        if idx <= gamma_ad {
            mx.trim_ad(bound, params.trim_thresh(max_score_in_ad, max_score));

            if bound.is_empty() {
                break;
            }
        }

        cloud.advance_reverse();
    }

    cloud.ad_start += 1;

    (1..=cell_start.seq_idx).rev().for_each(|seq_idx| {
        mx[(N, seq_idx)] = log_sum!(
            mx[(N, seq_idx + 1)] + prf[BackgroundLoop(N)],
            mx[(B, seq_idx)] + prf.special_transitions[Profile::N_IDX][Profile::SPECIAL_MOVE_IDX]
        );
    });
    CloudSearchResults {
        max_score: Nats(max_score),
        max_score_within: Nats(max_score_within),
        num_cells_computed: num_cells,
    }
}

#[inline]
pub fn compute_forward_cells<M>(
    prf: &Profile,
    seq: &Sequence,
    mx: &mut M,
    prf_idx: usize,
    seq_idx: usize,
) where
    M: NewDpMatrix,
{
    let residue = AminoAcid::from_u8(seq.digital_bytes[seq_idx]);

    // match state
    let m_dest = M(prf_idx);
    let m_cell = (m_dest, seq_idx);

    let m_src = M(prf_idx - 1);
    let i_src = I(prf_idx - 1);
    let d_src = D(prf_idx - 1);

    let m_src_cell = (m_src, seq_idx - 1);
    let i_src_cell = (i_src, seq_idx - 1);
    let d_src_cell = (d_src, seq_idx - 1);

    mx[m_cell] = log_sum!(
        mx[m_src_cell] + prf[CoreToCore(m_src, m_dest)],
        mx[i_src_cell] + prf[CoreToCore(i_src, m_dest)],
        mx[d_src_cell] + prf[CoreToCore(d_src, m_dest)],
        mx[(B, seq_idx - 1)] + prf.transition_score(Profile::B_M_IDX, prf_idx - 1)
    ) + prf[Emission(m_dest, residue)];

    // insert state
    let i_dest = I(prf_idx);
    let i_cell = (i_dest, seq_idx);

    let m_src = M(prf_idx);
    let i_src = I(prf_idx);

    let m_src_cell = (m_src, seq_idx - 1);
    let i_src_cell = (i_src, seq_idx - 1);

    mx[i_cell] = log_sum!(
        mx[m_src_cell] + prf[CoreToCore(m_src, i_dest)],
        mx[i_src_cell] + prf[CoreToCore(i_src, i_dest)]
    ) + prf[Emission(i_dest, residue)];

    // delete state
    let d_dest = D(prf_idx);
    let d_cell = (d_dest, seq_idx);

    let m_src = M(prf_idx - 1);
    let d_src = D(prf_idx - 1);

    let m_src_cell = (m_src, seq_idx);
    let d_src_cell = (d_src, seq_idx);

    mx[d_cell] = log_sum!(
        mx[m_src_cell] + prf[CoreToCore(m_src, d_dest)],
        mx[d_src_cell] + prf[CoreToCore(d_src, d_dest)]
    );

    mx[(E, seq_idx)] = log_sum!(mx[m_cell], mx[d_cell], mx[(E, seq_idx)])
}

pub fn cloud_search_forward<M>(
    prf: &Profile,
    seq: &Sequence,
    seed: &Seed,
    mx: &mut M,
    params: &CloudSearchParams,
    cloud: &mut Cloud,
) -> CloudSearchResults
where
    M: NewCloudMatrix,
{
    let mut num_cells = 0usize;
    let mut max_score = -f32::INFINITY;
    let mut max_score_within = -f32::INFINITY;

    let ad_start = seed.seq_start + seed.prf_start;
    let seed_ad_end = seed.seq_end + seed.prf_end;
    let gamma_ad = ad_start + params.gamma;
    let max_ad = seq.length + prf.length;

    let cell_start = Cell {
        prf_idx: seed.prf_start,
        seq_idx: seed.seq_start,
    };

    cloud.append(Bound(cell_start, cell_start));

    // precompute background (N) scores
    let seq_idx = cell_start.seq_idx - 1;
    let num_loops = seq_idx;
    mx[(N, seq_idx)] = num_loops as f32 * prf[BackgroundLoop(N)];
    mx[(B, seq_idx)] =
        mx[(N, seq_idx)] + prf.special_transitions[Profile::N_IDX][Profile::SPECIAL_MOVE_IDX];
    (cell_start.seq_idx..=seq.length).for_each(|seq_idx| {
        mx[(N, seq_idx)] = mx[(N, seq_idx - 1)] + prf[BackgroundLoop(N)];
        mx[(B, seq_idx)] =
            mx[(N, seq_idx)] + prf.special_transitions[Profile::N_IDX][Profile::SPECIAL_MOVE_IDX];
    });

    for idx in ad_start..=max_ad {
        let co_located_bound = &cloud[Ad(idx.saturating_sub(3))];
        mx.reset_ad(co_located_bound);

        let bound = &mut cloud[Ad(idx)];
        let mut max_score_in_ad = -f32::INFINITY;
        bound.iter().for_each(|c| {
            compute_forward_cells(prf, seq, mx, c.prf_idx, c.seq_idx);
            max_score_in_ad = max_f32!(
                max_score_in_ad,
                mx[c.m_cell()],
                mx[c.i_cell()],
                mx[c.d_cell()]
            );
            num_cells += 1;
        });

        max_score.max_assign(max_score_in_ad);
        if idx < seed_ad_end {
            max_score_within.max_assign(max_score)
        }

        if idx >= gamma_ad {
            mx.trim_ad(bound, params.trim_thresh(max_score_in_ad, max_score));

            if bound.is_empty() {
                break;
            }
        }

        cloud.advance_forward();
    }

    cloud.ad_end -= 1;

    (cell_start.seq_idx..=seq.length).for_each(|seq_idx| {
        mx[(C, seq_idx)] = log_sum!(
            mx[(C, seq_idx - 1)] + prf[BackgroundLoop(N)],
            mx[(E, seq_idx)] + prf.special_transitions[Profile::E_IDX][Profile::SPECIAL_MOVE_IDX]
        );
    });

    CloudSearchResults {
        max_score: Nats(max_score),
        max_score_within: Nats(max_score_within),
        num_cells_computed: num_cells,
    }
}
