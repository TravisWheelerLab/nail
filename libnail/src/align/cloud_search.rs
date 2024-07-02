use crate::align::structs::{AntiDiagonal, AntiDiagonalBounds, CloudMatrixLinear, Seed};
use crate::log_sum;
use crate::max_f32;
use crate::structs::{Profile, Sequence};
use crate::util::log_add;

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

pub struct CloudSearchScores {
    pub max_score: Nats,
    pub max_score_within: Nats,
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
    //   rotated linear orientation:
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
                + profile.transition_score(Profile::MATCH_TO_MATCH_IDX, profile_idx)
                + profile.match_score(previous_target_character, profile_idx + 1),
            cloud_matrix.get_insert(insert_source_row_idx, profile_idx)
                + profile.transition_score(Profile::MATCH_TO_INSERT_IDX, profile_idx)
                + profile.insert_score(previous_target_character, profile_idx),
            cloud_matrix.get_delete(delete_source_row_idx, profile_idx + 1)
                + profile.transition_score(Profile::MATCH_TO_DELETE_IDX, profile_idx)
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
                + profile.transition_score(Profile::INSERT_TO_MATCH_IDX, profile_idx)
                + profile.match_score(previous_target_character, profile_idx + 1),
            cloud_matrix.get_insert(insert_source_row_idx, profile_idx)
                + profile.transition_score(Profile::INSERT_TO_INSERT_IDX, profile_idx)
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
                + profile.transition_score(Profile::DELETE_TO_MATCH_IDX, profile_idx)
                + profile.match_score(previous_target_character, profile_idx + 1),
            cloud_matrix.get_delete(delete_source_row_idx, profile_idx + 1)
                + profile.transition_score(Profile::DELETE_TO_DELETE_IDX, profile_idx)
        ),
    );
}

pub fn cloud_search_backward(
    profile: &Profile,
    target: &Sequence,
    seed: &Seed,
    cloud_matrix: &mut CloudMatrixLinear,
    params: &CloudSearchParams,
    bounds: &mut AntiDiagonalBounds,
) -> CloudSearchScores {
    // the highest score we've seen overall
    let mut max_score = -f32::INFINITY;
    // the highest score we see before we pass the end seed point
    let mut max_score_within = -f32::INFINITY;

    let target_end = seed.target_end.min(target.length - 1);
    let profile_end = seed.profile_end.min(profile.length - 1);

    let first_anti_diagonal_idx = target_end + profile_end;
    let seed_start_anti_diagonal_idx = seed.target_start + seed.profile_start;
    let gamma_anti_diagonal_idx = first_anti_diagonal_idx - params.gamma;
    let min_anti_diagonal_idx = 0usize;

    let first_cloud_matrix_row_idx = first_anti_diagonal_idx % 3;

    // setting the scores to 0 is like setting
    // the log odds ratio to 1 since log(0) = 1
    cloud_matrix.set_match(first_cloud_matrix_row_idx, profile_end, 0.0);
    cloud_matrix.set_insert(first_cloud_matrix_row_idx, profile_end, 0.0);
    cloud_matrix.set_delete(first_cloud_matrix_row_idx, profile_end, 0.0);

    // the first bound is the ending position
    bounds.max_anti_diagonal_idx = first_anti_diagonal_idx;
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
                bounds.min_anti_diagonal_idx = anti_diagonal_idx + 1;
                break;
            }
            PruneStatus::PartiallyPruned => {
                bounds.min_anti_diagonal_idx = anti_diagonal_idx;
                continue;
            }
        }
    }

    CloudSearchScores {
        max_score: Nats(max_score),
        max_score_within: Nats(max_score_within),
    }
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
    //   rotated linear orientation:
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
                + profile.transition_score(Profile::MATCH_TO_MATCH_IDX, source_profile_idx),
            cloud_matrix.get_insert(source_row_idx, source_profile_idx)
                + profile.transition_score(Profile::INSERT_TO_MATCH_IDX, source_profile_idx),
            cloud_matrix.get_delete(source_row_idx, source_profile_idx)
                + profile.transition_score(Profile::DELETE_TO_MATCH_IDX, source_profile_idx)
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
                + profile.transition_score(Profile::MATCH_TO_INSERT_IDX, profile_idx),
            cloud_matrix.get_insert(source_row_idx, profile_idx)
                + profile.transition_score(Profile::INSERT_TO_INSERT_IDX, profile_idx)
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
                + profile.transition_score(Profile::MATCH_TO_DELETE_IDX, source_profile_idx),
            cloud_matrix.get_delete(source_row_idx, source_profile_idx)
                + profile.transition_score(Profile::DELETE_TO_DELETE_IDX, source_profile_idx)
        ),
    );
}

pub fn cloud_search_forward(
    profile: &Profile,
    target: &Sequence,
    seed: &Seed,
    cloud_matrix: &mut CloudMatrixLinear,
    params: &CloudSearchParams,
    bounds: &mut AntiDiagonalBounds,
) -> CloudSearchScores {
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

    let first_anti_diagonal_idx = seed.target_start + seed.profile_start;
    let seed_end_anti_diagonal_idx = seed.target_end + seed.profile_end;
    let gamma_anti_diagonal_idx = first_anti_diagonal_idx + params.gamma;
    let max_anti_diagonal_idx = target.length + profile.length;

    let first_cloud_matrix_row_idx = first_anti_diagonal_idx % 3;
    // setting the scores to 0 is like setting
    // the log odds ratio to 1, since log(0) = 1
    cloud_matrix.set_match(first_cloud_matrix_row_idx, seed.profile_start, 0.0);
    cloud_matrix.set_insert(first_cloud_matrix_row_idx, seed.profile_start, 0.0);
    cloud_matrix.set_delete(first_cloud_matrix_row_idx, seed.profile_start, 0.0);

    // the first bound is just the starting cell
    bounds.min_anti_diagonal_idx = first_anti_diagonal_idx;
    bounds.set(
        first_anti_diagonal_idx,
        seed.target_start,
        seed.profile_start,
        seed.target_start,
        seed.profile_start,
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
                bounds.max_anti_diagonal_idx = anti_diagonal_idx - 1;
                break;
            }
            PruneStatus::PartiallyPruned => {
                bounds.max_anti_diagonal_idx = anti_diagonal_idx;
                continue;
            }
        }
    }

    CloudSearchScores {
        max_score: Nats(max_score),
        max_score_within: Nats(max_score_within),
    }
}
