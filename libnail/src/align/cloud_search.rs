use crate::align::structs::{AntiDiagonal, Cloud, CloudMatrixLinear, Seed};
use crate::log_sum;
use crate::max_f32;
use crate::structs::profile::{AminoAcid, BackgroundLoop, CoreToCore, Emission};
use crate::structs::{Profile, Sequence};
use crate::util::log_add;

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

pub enum PruneStatus {
    FullyPruned,
    PartiallyPruned,
}

#[derive(Default, Debug, PartialEq)]
pub struct CloudSearchResults {
    pub max_score: Nats,
    pub max_score_within: Nats,
    pub num_cells_computed: usize,
}

#[inline]
pub fn prune_and_scrub(
    bound: &mut AntiDiagonal,
    cloud_matrix: &mut CloudMatrixLinear,
    row_idx: usize,
    alpha: f32,
    beta: f32,
    overall_max: &mut f32,
) -> PruneStatus {
    let mut current_max = -f32::INFINITY;

    // TODO: would we get better performance if we keep track of the max values here?
    for profile_idx in bound.left_profile_idx..=bound.right_profile_idx {
        let max_score = max_f32!(
            cloud_matrix.get_match(row_idx, profile_idx),
            cloud_matrix.get_insert(row_idx, profile_idx),
            cloud_matrix.get_delete(row_idx, profile_idx)
        );
        current_max = current_max.max(max_score);
        *overall_max = overall_max.max(max_score);
    }

    let alpha_thresh = current_max - alpha;
    let beta_thresh = *overall_max - beta;

    for profile_idx in bound.left_profile_idx..=bound.right_profile_idx {
        let max_score = max_f32!(
            cloud_matrix.get_match(row_idx, profile_idx),
            cloud_matrix.get_insert(row_idx, profile_idx),
            cloud_matrix.get_delete(row_idx, profile_idx)
        );
        if max_score > alpha_thresh && max_score > beta_thresh {
            break;
        } else {
            bound.left_profile_idx += 1;
            bound.left_target_idx -= 1;
            cloud_matrix.set_match(row_idx, profile_idx, -f32::INFINITY);
            cloud_matrix.set_insert(row_idx, profile_idx, -f32::INFINITY);
            cloud_matrix.set_delete(row_idx, profile_idx, -f32::INFINITY);
        }
    }

    for profile_idx in (bound.left_profile_idx..=bound.right_profile_idx).rev() {
        let max_score = max_f32!(
            cloud_matrix.get_match(row_idx, profile_idx),
            cloud_matrix.get_insert(row_idx, profile_idx),
            cloud_matrix.get_delete(row_idx, profile_idx)
        );
        if max_score > alpha_thresh && max_score > beta_thresh {
            break;
        } else {
            bound.right_profile_idx -= 1;
            bound.right_target_idx += 1;
            cloud_matrix.set_match(row_idx, profile_idx, -f32::INFINITY);
            cloud_matrix.set_insert(row_idx, profile_idx, -f32::INFINITY);
            cloud_matrix.set_delete(row_idx, profile_idx, -f32::INFINITY);
        }
    }

    if bound.was_pruned() {
        bound.reset();
        PruneStatus::FullyPruned
    } else {
        PruneStatus::PartiallyPruned
    }
}

#[inline]
pub fn scrub_co_located(
    current_bound: &AntiDiagonal,
    co_located_bound: &AntiDiagonal,
    cloud_matrix: &mut CloudMatrixLinear,
    cloud_matrix_row_idx: usize,
) {
    let left_scrub_amount = current_bound
        .left_profile_idx
        .saturating_sub(co_located_bound.left_profile_idx);

    let right_scrub_amount = co_located_bound
        .right_profile_idx
        .saturating_sub(current_bound.right_profile_idx);

    let left_scrub_start = current_bound.left_profile_idx - left_scrub_amount;
    let right_scrub_end = current_bound.right_profile_idx + right_scrub_amount;

    // TODO: double check these bounds
    for profile_idx in left_scrub_start..current_bound.left_profile_idx {
        cloud_matrix.set_match(cloud_matrix_row_idx, profile_idx, -f32::INFINITY);
        cloud_matrix.set_insert(cloud_matrix_row_idx, profile_idx, -f32::INFINITY);
        cloud_matrix.set_delete(cloud_matrix_row_idx, profile_idx, -f32::INFINITY);
    }

    for profile_idx in (current_bound.right_profile_idx + 1)..=right_scrub_end {
        cloud_matrix.set_match(cloud_matrix_row_idx, profile_idx, -f32::INFINITY);
        cloud_matrix.set_insert(cloud_matrix_row_idx, profile_idx, -f32::INFINITY);
        cloud_matrix.set_delete(cloud_matrix_row_idx, profile_idx, -f32::INFINITY);
    }
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

#[inline]
pub fn compute_backward_cell(
    target: &Sequence,
    profile: &Profile,
    cloud_matrix: &mut CloudMatrixLinear,
    cloud_matrix_row_idx: usize,
    target_idx: usize,
    profile_idx: usize,
) {
    let previous_target_character = target.digital_bytes[target_idx + 1] as usize;

    //   *: the cell we are computing                     (target_idx    , profile_idx    )
    //   M: the cell of the source match state component  (target_idx + 1, profile_idx + 1)
    //   D: the cell of the source delete state component (target_idx    , profile_idx + 1)
    //   I: the cell of the source match state component  (target_idx + 1, profile_idx    )
    //
    //   classic orientation:
    //       p r o f i l e
    //     t - - - - - - - -
    //     a - - - - - - - -
    //     r - - - * D - - -
    //     g - - - I M - - -
    //     e - - - - - - - -
    //     t - - - - - - - -
    //
    //   linear orientation:
    //       p r o f i l e
    //     r - - - * - - - -
    //     o - - - I D - - -
    //     w - - - - M - - -
    //
    let match_source_row_idx = (cloud_matrix_row_idx + 2) % 3;
    let insert_source_row_idx = (cloud_matrix_row_idx + 1) % 3;
    let delete_source_row_idx = (cloud_matrix_row_idx + 1) % 3;

    // match state
    //
    //  * D    >   * -
    //  I M    >   I D
    //         >   - M
    //
    cloud_matrix.set_match(
        cloud_matrix_row_idx,
        profile_idx,
        log_sum!(
            cloud_matrix.get_match(match_source_row_idx, profile_idx + 1)
                + profile.transition_score(Profile::M_M_IDX, profile_idx)
                + profile.match_score(previous_target_character, profile_idx + 1),
            cloud_matrix.get_insert(insert_source_row_idx, profile_idx)
                + profile.transition_score(Profile::M_I_IDX, profile_idx)
                + profile.insert_score(previous_target_character, profile_idx),
            cloud_matrix.get_delete(delete_source_row_idx, profile_idx + 1)
                + profile.transition_score(Profile::M_D_IDX, profile_idx)
        ),
    );

    // insert state
    //
    //  * -   >   * -
    //  I M   >   I -
    //        >   - M
    //
    cloud_matrix.set_insert(
        cloud_matrix_row_idx,
        profile_idx,
        log_sum!(
            cloud_matrix.get_match(match_source_row_idx, profile_idx + 1)
                + profile.transition_score(Profile::I_M_IDX, profile_idx)
                + profile.match_score(previous_target_character, profile_idx + 1),
            cloud_matrix.get_insert(insert_source_row_idx, profile_idx)
                + profile.transition_score(Profile::I_I_IDX, profile_idx)
                + profile.insert_score(previous_target_character, profile_idx)
        ),
    );

    // delete state
    //
    //   * D    >    * -
    //   - M    >    - D
    //          >    - M
    //
    cloud_matrix.set_delete(
        cloud_matrix_row_idx,
        profile_idx,
        log_sum!(
            cloud_matrix.get_match(match_source_row_idx, profile_idx + 1)
                + profile.transition_score(Profile::D_M_IDX, profile_idx)
                + profile.match_score(previous_target_character, profile_idx + 1),
            cloud_matrix.get_delete(delete_source_row_idx, profile_idx + 1)
                + profile.transition_score(Profile::D_D_IDX, profile_idx)
        ),
    );
}

pub fn cloud_search_backward2<M>(
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

        max_score = max_score.max(max_score_in_ad);
        if idx > seed_ad_start {
            max_score_within = max_score_within.max(max_score)
        }

        if idx <= gamma_ad {
            let trim_thresh = (max_score_in_ad - params.alpha).max(max_score - params.beta);
            let trim_fn = |mx: &mut M, c: &Cell| {
                let max = max_f32!(mx[c.m_cell()], mx[c.i_cell()], mx[c.d_cell()]);
                if max < trim_thresh {
                    mx[c.m_cell()] = -f32::INFINITY;
                    mx[c.i_cell()] = -f32::INFINITY;
                    mx[c.d_cell()] = -f32::INFINITY;
                    true
                } else {
                    false
                }
            };

            let left_trim = bound.iter().take_while(|c| trim_fn(mx, c)).count();
            bound.0.prf_idx -= left_trim;
            bound.0.seq_idx += left_trim;

            let right_trim = bound.iter().rev().take_while(|c| trim_fn(mx, c)).count();
            bound.1.prf_idx += right_trim;
            bound.1.seq_idx -= right_trim;

            if bound.is_empty() {
                break;
            }
        }

        cloud.advance_reverse();
    }

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

pub fn cloud_search_backward(
    profile: &Profile,
    target: &Sequence,
    seed: &Seed,
    cloud_matrix: &mut CloudMatrixLinear,
    params: &CloudSearchParams,
    bounds: &mut Cloud,
) -> CloudSearchResults {
    let mut num_cells_computed = 0;
    // the highest score we've seen overall
    let mut max_score = -f32::INFINITY;
    // the highest score we see before we pass the end seed point
    let mut max_score_within = -f32::INFINITY;

    let target_end = seed.seq_end.min(target.length - 1);
    let profile_end = seed.prf_end.min(profile.length - 1);

    let first_anti_diagonal_idx = target_end + profile_end;
    let seed_start_anti_diagonal_idx = seed.seq_start + seed.prf_start;
    let gamma_anti_diagonal_idx = first_anti_diagonal_idx - params.gamma;
    let min_anti_diagonal_idx = 0usize;

    let first_cloud_matrix_row_idx = first_anti_diagonal_idx % 3;

    // setting the scores to 0 is like setting
    // the log odds ratio to 1 since log(0) = 1
    cloud_matrix.set_match(first_cloud_matrix_row_idx, profile_end, 0.0);
    cloud_matrix.set_insert(first_cloud_matrix_row_idx, profile_end, 0.0);
    cloud_matrix.set_delete(first_cloud_matrix_row_idx, profile_end, 0.0);

    // the first bound is the ending position
    bounds.ad_end = first_anti_diagonal_idx;
    bounds.set(
        first_anti_diagonal_idx,
        target_end,
        profile_end,
        target_end,
        profile_end,
    );

    for anti_diagonal_idx in (gamma_anti_diagonal_idx..first_anti_diagonal_idx).rev() {
        let previous_bound = bounds.get(anti_diagonal_idx + 1);

        bounds.set(
            anti_diagonal_idx,
            previous_bound.left_target_idx,
            previous_bound.left_profile_idx,
            previous_bound.right_target_idx,
            previous_bound.right_profile_idx,
        );

        let current_bound = bounds.get_mut(anti_diagonal_idx);
        // edge guards
        if current_bound.left_profile_idx > 1 {
            // if we haven't hit the start of the target (top row),
            // we'll plan to move the left bound up
            current_bound.left_profile_idx -= 1;
        } else {
            // otherwise we'll move to the left
            current_bound.left_target_idx -= 1;
        }

        if current_bound.right_target_idx > 1 {
            // if we haven't hit the start of the profile (first column),
            // we'll plan to move the right bound to the left
            current_bound.right_target_idx -= 1;
        } else {
            // otherwise we'll move up
            current_bound.right_profile_idx -= 1;
        }

        let cloud_matrix_row_idx = anti_diagonal_idx % 3;

        for (target_idx, profile_idx) in current_bound.cell_zip() {
            num_cells_computed += 1;
            compute_backward_cell(
                target,
                profile,
                cloud_matrix,
                cloud_matrix_row_idx,
                target_idx,
                profile_idx,
            );
        }
    }

    // main recursion:
    //  - compute the next full anti-diagonal
    //  - prune with max-in-anti-diagonal rule: discard any value < max_in_current_diagonal - alpha (default: alpha = 12)
    //  - prune with x-drop rule: discard anything value < max_in_all_diagonals - beta (default: beta = 20)
    for anti_diagonal_idx in (min_anti_diagonal_idx..gamma_anti_diagonal_idx).rev() {
        let previous_bound = bounds.get(anti_diagonal_idx + 1);
        let three_forward_bound = bounds.get(anti_diagonal_idx + 3).clone();

        bounds.set(
            anti_diagonal_idx,
            previous_bound.left_target_idx,
            previous_bound.left_profile_idx,
            previous_bound.right_target_idx,
            previous_bound.right_profile_idx,
        );

        let current_bound = bounds.get_mut(anti_diagonal_idx);
        // edge guards
        if current_bound.left_profile_idx > 1 {
            // if we haven't hit the start of the target (top row),
            // we'll plan to move the left bound up
            current_bound.left_profile_idx -= 1;
        } else {
            // otherwise we'll move to the left
            current_bound.left_target_idx -= 1;
        }

        if current_bound.right_target_idx > 1 {
            // if we haven't hit the start of the profile (first column),
            // we'll plan to move the right bound to the left
            current_bound.right_target_idx -= 1;
        } else {
            // otherwise we'll move up
            current_bound.right_profile_idx -= 1;
        }

        let cloud_matrix_row_idx = anti_diagonal_idx % 3;

        // 3-forward scrub
        scrub_co_located(
            current_bound,
            &three_forward_bound,
            cloud_matrix,
            cloud_matrix_row_idx,
        );

        for (target_idx, profile_idx) in current_bound.cell_zip() {
            num_cells_computed += 1;
            compute_backward_cell(
                target,
                profile,
                cloud_matrix,
                cloud_matrix_row_idx,
                target_idx,
                profile_idx,
            );
        }

        let prune_status = prune_and_scrub(
            current_bound,
            cloud_matrix,
            cloud_matrix_row_idx,
            params.alpha,
            params.beta,
            &mut max_score,
        );

        if anti_diagonal_idx >= seed_start_anti_diagonal_idx {
            max_score_within = max_score_within.max(max_score);
        }

        match prune_status {
            PruneStatus::FullyPruned => {
                bounds.ad_start = anti_diagonal_idx + 1;
                break;
            }
            PruneStatus::PartiallyPruned => {
                bounds.ad_start = anti_diagonal_idx;
                continue;
            }
        }
    }

    CloudSearchResults {
        max_score: Nats(max_score),
        max_score_within: Nats(max_score_within),
        num_cells_computed,
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

#[inline]
pub fn compute_forward_cell(
    target: &Sequence,
    profile: &Profile,
    cloud_matrix: &mut CloudMatrixLinear,
    cloud_matrix_row_idx: usize,
    target_idx: usize,
    profile_idx: usize,
) {
    let current_target_character = target.digital_bytes[target_idx];

    // match state
    // note: begin to match excluded here
    //   *: the cell we are computing (target_idx    , profile_idx    )
    //   S: the source cell           (target_idx - 1, profile_idx - 1)
    //
    //   classic orientation:
    //       p r o f i l e
    //     t - - - - - - - -
    //     a - - - - - - - -
    //     r - - - S - - - -
    //     g - - - - * - - -
    //     e - - - - - - - -
    //     t - - - - - - - -
    //
    //   linear orientation:
    //       p r o f i l e
    //     r - - - S - - - -
    //     o - - - - - - - -
    //     w - - - - * - - -
    //
    let mut source_row_idx = (cloud_matrix_row_idx + 1) % 3;
    let mut source_profile_idx = profile_idx - 1;
    cloud_matrix.set_match(
        cloud_matrix_row_idx,
        profile_idx,
        log_sum!(
            cloud_matrix.get_match(source_row_idx, source_profile_idx)
                + profile.transition_score(Profile::M_M_IDX, source_profile_idx),
            cloud_matrix.get_insert(source_row_idx, source_profile_idx)
                + profile.transition_score(Profile::I_M_IDX, source_profile_idx),
            cloud_matrix.get_delete(source_row_idx, source_profile_idx)
                + profile.transition_score(Profile::D_M_IDX, source_profile_idx)
        ) + profile.match_score(current_target_character as usize, profile_idx),
    );

    // insert state
    //
    //  - S    >   - -
    //  - *    >   - S
    //         >   - *
    //

    source_row_idx = (cloud_matrix_row_idx + 2) % 3;
    cloud_matrix.set_insert(
        cloud_matrix_row_idx,
        profile_idx,
        log_sum!(
            cloud_matrix.get_match(source_row_idx, profile_idx)
                + profile.transition_score(Profile::M_I_IDX, profile_idx),
            cloud_matrix.get_insert(source_row_idx, profile_idx)
                + profile.transition_score(Profile::I_I_IDX, profile_idx)
        ) + profile.insert_score(current_target_character as usize, profile_idx),
    );

    // delete state
    //
    //  - -    >   - -
    //  S *    >   S -
    //         >   - *
    //
    source_row_idx = (cloud_matrix_row_idx + 2) % 3;
    source_profile_idx = profile_idx - 1;
    cloud_matrix.set_delete(
        cloud_matrix_row_idx,
        profile_idx,
        log_sum!(
            cloud_matrix.get_match(source_row_idx, source_profile_idx)
                + profile.transition_score(Profile::M_D_IDX, source_profile_idx),
            cloud_matrix.get_delete(source_row_idx, source_profile_idx)
                + profile.transition_score(Profile::D_D_IDX, source_profile_idx)
        ),
    );
}

pub fn cloud_search_forward(
    profile: &Profile,
    target: &Sequence,
    seed: &Seed,
    cloud_matrix: &mut CloudMatrixLinear,
    params: &CloudSearchParams,
    bounds: &mut Cloud,
) -> CloudSearchResults {
    let mut num_cells_computed = 0;
    // the highest score we've seen overall
    let mut max_score = -f32::INFINITY;
    // the highest score we see before we pass the end seed point
    let mut max_score_within = -f32::INFINITY;

    // the first valid anti_diagonal_idx is 2
    //
    //                    x x x x x x x x x x
    //                    x 2 - - - - - - - -
    // target_start = 2 > x - - 5 6 7 8 9 - -
    //                    x - - 6 7 8 9 - - -
    //                    x - - 7 8 9 - - - -
    //                    x - - 8 9 - - - - -
    //                    x - - 9 - - - - - -
    //                    x - - - - - - - - -
    //                          ^
    //                   profile_start = 3

    let first_anti_diagonal_idx = seed.seq_start + seed.prf_start;
    let seed_end_anti_diagonal_idx = seed.seq_end + seed.prf_end;
    let gamma_anti_diagonal_idx = first_anti_diagonal_idx + params.gamma;
    let max_anti_diagonal_idx = target.length + profile.length;

    let first_cloud_matrix_row_idx = first_anti_diagonal_idx % 3;
    // setting the scores to 0 is like setting
    // the log odds ratio to 1, since log(0) = 1
    cloud_matrix.set_match(first_cloud_matrix_row_idx, seed.prf_start, 0.0);
    cloud_matrix.set_insert(first_cloud_matrix_row_idx, seed.prf_start, 0.0);
    cloud_matrix.set_delete(first_cloud_matrix_row_idx, seed.prf_start, 0.0);

    // the first bound is just the starting cell
    bounds.ad_start = first_anti_diagonal_idx;
    bounds.set(
        first_anti_diagonal_idx,
        seed.seq_start,
        seed.prf_start,
        seed.seq_start,
        seed.prf_start,
    );

    for anti_diagonal_idx in (first_anti_diagonal_idx + 1)..gamma_anti_diagonal_idx {
        let previous_bound = bounds.get(anti_diagonal_idx - 1);

        bounds.set(
            anti_diagonal_idx,
            previous_bound.left_target_idx,
            previous_bound.left_profile_idx,
            previous_bound.right_target_idx,
            previous_bound.right_profile_idx,
        );

        let current_bound = bounds.get_mut(anti_diagonal_idx);

        // edge guards
        if current_bound.left_target_idx < target.length {
            // if we haven't hit the end of the target (bottom row),
            // we'll plan to move the left bound down 1
            current_bound.left_target_idx += 1;
        } else {
            // otherwise we'll move to the right 1
            current_bound.left_profile_idx += 1;
        }

        if current_bound.right_profile_idx < profile.length {
            // if we haven't hit the end of the profile (last column),
            // we'll plan to move the right bound to the right 1
            current_bound.right_profile_idx += 1;
        } else {
            // otherwise we'll move down 1
            current_bound.right_target_idx += 1;
        }

        let cloud_matrix_row_idx = anti_diagonal_idx % 3;

        for (target_idx, profile_idx) in current_bound.cell_zip() {
            num_cells_computed += 1;
            compute_forward_cell(
                target,
                profile,
                cloud_matrix,
                cloud_matrix_row_idx,
                target_idx,
                profile_idx,
            );
        }
    }
    // main recursion:
    //  - compute the next full anti-diagonal
    //  - prune with max-in-anti-diagonal rule: discard any value < max_in_current_diagonal - alpha (default: alpha = 12)
    //  - prune with x-drop rule: discard anything value < max_in_all_diagonals - beta (default: beta = 20)
    for anti_diagonal_idx in gamma_anti_diagonal_idx..=max_anti_diagonal_idx {
        let previous_bound = bounds.get(anti_diagonal_idx - 1);
        let three_back_bound = bounds.get(anti_diagonal_idx - 3).clone();

        bounds.set(
            anti_diagonal_idx,
            previous_bound.left_target_idx,
            previous_bound.left_profile_idx,
            previous_bound.right_target_idx,
            previous_bound.right_profile_idx,
        );

        let current_bound = bounds.get_mut(anti_diagonal_idx);

        // edge guards
        if current_bound.left_target_idx < target.length {
            // if we haven't hit the end of the target (bottom row),
            // we'll plan to move the left bound down 1
            current_bound.left_target_idx += 1;
        } else {
            // otherwise we'll move to the right 1
            current_bound.left_profile_idx += 1;
        }

        if current_bound.right_profile_idx < profile.length {
            // if we haven't hit the end of the profile (last column),
            // we'll plan to move the right bound to the right 1
            current_bound.right_profile_idx += 1;
        } else {
            // otherwise we'll move down 1
            current_bound.right_target_idx += 1;
        }

        let cloud_matrix_row_idx = anti_diagonal_idx % 3;

        // 3-back scrub
        scrub_co_located(
            current_bound,
            &three_back_bound,
            cloud_matrix,
            cloud_matrix_row_idx,
        );

        for (target_idx, profile_idx) in current_bound.cell_zip() {
            num_cells_computed += 1;
            compute_forward_cell(
                target,
                profile,
                cloud_matrix,
                cloud_matrix_row_idx,
                target_idx,
                profile_idx,
            );
        }

        let prune_status = prune_and_scrub(
            current_bound,
            cloud_matrix,
            cloud_matrix_row_idx,
            params.alpha,
            params.beta,
            &mut max_score,
        );

        if anti_diagonal_idx <= seed_end_anti_diagonal_idx {
            max_score_within = max_score_within.max(max_score);
        }

        match prune_status {
            PruneStatus::FullyPruned => {
                bounds.ad_end = anti_diagonal_idx - 1;
                break;
            }
            PruneStatus::PartiallyPruned => {
                bounds.ad_end = anti_diagonal_idx;
                continue;
            }
        }
    }

    CloudSearchResults {
        max_score: Nats(max_score),
        max_score_within: Nats(max_score_within),
        num_cells_computed,
    }
}

pub fn cloud_search_forward2<M>(
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

        max_score = max_score.max(max_score_in_ad);
        if idx < seed_ad_end {
            max_score_within = max_score_within.max(max_score)
        }

        if idx >= gamma_ad {
            let trim_thresh = (max_score_in_ad - params.alpha).max(max_score - params.beta);
            let trim_fn = |mx: &mut M, c: &Cell| {
                let max = max_f32!(mx[c.m_cell()], mx[c.i_cell()], mx[c.d_cell()]);
                if max < trim_thresh {
                    mx[c.m_cell()] = -f32::INFINITY;
                    mx[c.i_cell()] = -f32::INFINITY;
                    mx[c.d_cell()] = -f32::INFINITY;
                    true
                } else {
                    false
                }
            };

            let left_trim = bound.iter().take_while(|c| trim_fn(mx, c)).count();
            bound.0.prf_idx -= left_trim;
            bound.0.seq_idx += left_trim;

            let right_trim = bound.iter().rev().take_while(|c| trim_fn(mx, c)).count();
            bound.1.prf_idx += right_trim;
            bound.1.seq_idx -= right_trim;

            if bound.is_empty() {
                break;
            }
        }

        cloud.advance_forward();
    }

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
