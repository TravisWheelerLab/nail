use crate::align::bounded::cloud_search_common::PruneStatus;
use crate::align::bounded::structs::{CloudBoundGroup, CloudMatrixLinear, CloudSearchParams, Seed};
use crate::align::bounded::{prune_and_scrub, scrub_co_located};
use crate::log_sum;
use crate::structs::{Profile, Sequence};
use crate::util::log_add;

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
    bounds: &mut CloudBoundGroup,
) {
    // the highest score we've seen overall
    let mut overall_max_score = -f32::INFINITY;

    let target_end = seed.target_end.min(target.length - 1);
    let profile_end = seed.profile_end.min(profile.length - 1);

    let first_anti_diagonal_idx = target_end + profile_end;
    let gamma_anti_diagonal_idx = first_anti_diagonal_idx - params.gamma;
    let min_anti_diagonal_idx = 0usize;

    let first_cloud_matrix_row_idx = first_anti_diagonal_idx % 3;

    // setting the scores to 0 is like setting
    // the log odds ratio to 1 since log(0) = 1
    cloud_matrix.set_match(first_cloud_matrix_row_idx, profile_end, 0.0);
    cloud_matrix.set_insert(first_cloud_matrix_row_idx, profile_end, 0.0);
    cloud_matrix.set_delete(first_cloud_matrix_row_idx, profile_end, 0.0);

    // the first bound is the ending position
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
            previous_bound.left_profile_idx - 1,
            previous_bound.right_target_idx - 1,
            previous_bound.right_profile_idx,
        );

        let current_bound = bounds.get(anti_diagonal_idx);
        let cloud_matrix_row_idx = anti_diagonal_idx % 3;

        for (target_idx, profile_idx) in current_bound.anti_diagonal_cell_zip() {
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

        for (target_idx, profile_idx) in current_bound.anti_diagonal_cell_zip() {
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
            &mut overall_max_score,
        );

        match prune_status {
            PruneStatus::FullyPruned => {
                bounds.min_anti_diagonal_idx += 1;
                break;
            }
            PruneStatus::PartiallyPruned => {
                continue;
            }
        }
    }
}
