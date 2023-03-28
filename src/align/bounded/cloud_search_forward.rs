use crate::align::bounded::cloud_search_common::{prune_and_scrub, scrub_co_located, PruneStatus};
use crate::align::bounded::structs::{CloudBoundGroup, CloudMatrixLinear, CloudSearchParams};
use crate::log_sum;
use crate::structs::{Profile, Sequence};
use crate::timing::time;
use crate::util::log_add;

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

#[funci::timed(timer = time)]
pub fn cloud_search_forward(
    profile: &Profile,
    target: &Sequence,
    cloud_matrix: &mut CloudMatrixLinear,
    params: &CloudSearchParams,
    bounds: &mut CloudBoundGroup,
) -> anyhow::Result<()> {
    // the highest score we've seen overall
    let mut overall_max_score = -f32::INFINITY;

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

    let first_anti_diagonal_idx = params.target_start + params.profile_start;
    let gamma_anti_diagonal_idx = first_anti_diagonal_idx + params.gamma;
    let max_anti_diagonal_idx = target.length + profile.length;

    let first_cloud_matrix_row_idx = first_anti_diagonal_idx % 3;
    // setting the scores to 0 is like setting
    // the log odds ratio to 1, since log(0) = 1
    cloud_matrix.set_match(first_cloud_matrix_row_idx, params.profile_start, 0.0);
    cloud_matrix.set_insert(first_cloud_matrix_row_idx, params.profile_start, 0.0);
    cloud_matrix.set_delete(first_cloud_matrix_row_idx, params.profile_start, 0.0);

    // the first bound is just the starting cell
    bounds.set(
        first_anti_diagonal_idx,
        params.target_start,
        params.profile_start,
        params.target_start,
        params.profile_start,
    );

    for anti_diagonal_idx in (first_anti_diagonal_idx + 1)..gamma_anti_diagonal_idx {
        let previous_bound = bounds.get(anti_diagonal_idx - 1);

        bounds.set(
            anti_diagonal_idx,
            previous_bound.left_target_idx + 1,
            previous_bound.left_profile_idx,
            previous_bound.right_target_idx,
            previous_bound.right_profile_idx + 1,
        );

        let current_bound = bounds.get(anti_diagonal_idx);
        let cloud_matrix_row_idx = anti_diagonal_idx % 3;

        for (target_idx, profile_idx) in current_bound.anti_diagonal_cell_zip() {
            compute_forward_cell(
                target,
                profile,
                cloud_matrix,
                cloud_matrix_row_idx,
                target_idx,
                profile_idx,
                // &mut debug,
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

        for (target_idx, profile_idx) in current_bound.anti_diagonal_cell_zip() {
            compute_forward_cell(
                target,
                profile,
                cloud_matrix,
                cloud_matrix_row_idx,
                target_idx,
                profile_idx,
                // &mut debug,
            );
        }

        let prune_status = prune_and_scrub(
            current_bound,
            cloud_matrix,
            cloud_matrix_row_idx,
            params.alpha,
            params.beta,
            &mut overall_max_score,
            // &mut debug,
        );

        match prune_status {
            PruneStatus::FullyPruned => {
                bounds.max_anti_diagonal_idx -= 1;
                break;
            }
            PruneStatus::PartiallyPruned => {
                continue;
            }
        }
    }

    Ok(())
}
