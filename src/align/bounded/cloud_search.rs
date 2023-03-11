use crate::align::bounded::structs::bound::{join_bounds, CloudBoundGroup};
use crate::align::bounded::structs::{CloudBound, CloudMatrixLinear, CloudSearchParams};
use crate::structs::profile::constants::{
    PROFILE_DELETE_TO_DELETE, PROFILE_DELETE_TO_MATCH, PROFILE_INSERT_TO_INSERT,
    PROFILE_INSERT_TO_MATCH, PROFILE_MATCH_TO_DELETE, PROFILE_MATCH_TO_INSERT,
    PROFILE_MATCH_TO_MATCH,
};
use crate::structs::{DpMatrix, Profile, Sequence};
use crate::util::{log_add, PrintMe};
use crate::viz::SodaJson;
use crate::{log_sum, max_f32};
use anyhow::Result;
use std::fs::File;
use std::io::BufWriter;

#[inline]
pub fn compute_forward_cell(
    target: &Sequence,
    profile: &Profile,
    cloud_matrix: &mut CloudMatrixLinear,
    cloud_matrix_row_idx: usize,
    target_idx: usize,
    profile_idx: usize,
    debug: &mut DpMatrix,
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
                + profile.transition_score(PROFILE_MATCH_TO_MATCH, source_profile_idx),
            cloud_matrix.get_insert(source_row_idx, source_profile_idx)
                + profile.transition_score(PROFILE_INSERT_TO_MATCH, source_profile_idx),
            cloud_matrix.get_delete(source_row_idx, source_profile_idx)
                + profile.transition_score(PROFILE_DELETE_TO_MATCH, source_profile_idx)
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
                + profile.transition_score(PROFILE_MATCH_TO_INSERT, profile_idx),
            cloud_matrix.get_insert(source_row_idx, profile_idx)
                + profile.transition_score(PROFILE_INSERT_TO_INSERT, profile_idx)
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
                + profile.transition_score(PROFILE_MATCH_TO_DELETE, source_profile_idx),
            cloud_matrix.get_delete(source_row_idx, source_profile_idx)
                + profile.transition_score(PROFILE_DELETE_TO_DELETE, source_profile_idx)
        ),
    );

    debug.set_match(
        target_idx,
        profile_idx,
        cloud_matrix.get_match(cloud_matrix_row_idx, profile_idx),
    );
    debug.set_insert(
        target_idx,
        profile_idx,
        cloud_matrix.get_insert(cloud_matrix_row_idx, profile_idx),
    );
    debug.set_delete(
        target_idx,
        profile_idx,
        cloud_matrix.get_delete(cloud_matrix_row_idx, profile_idx),
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
    debug: &mut DpMatrix,
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
                + profile.transition_score(PROFILE_MATCH_TO_MATCH, profile_idx)
                + profile.match_score(previous_target_character, profile_idx + 1),
            cloud_matrix.get_insert(insert_source_row_idx, profile_idx)
                + profile.transition_score(PROFILE_MATCH_TO_INSERT, profile_idx)
                + profile.insert_score(previous_target_character, profile_idx),
            cloud_matrix.get_delete(delete_source_row_idx, profile_idx + 1)
                + profile.transition_score(PROFILE_MATCH_TO_DELETE, profile_idx)
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
                + profile.transition_score(PROFILE_INSERT_TO_MATCH, profile_idx)
                + profile.match_score(previous_target_character, profile_idx + 1),
            cloud_matrix.get_insert(insert_source_row_idx, profile_idx)
                + profile.transition_score(PROFILE_INSERT_TO_INSERT, profile_idx)
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
                + profile.transition_score(PROFILE_DELETE_TO_MATCH, profile_idx)
                + profile.match_score(previous_target_character, profile_idx + 1),
            cloud_matrix.get_delete(delete_source_row_idx, profile_idx + 1)
                + profile.transition_score(PROFILE_DELETE_TO_DELETE, profile_idx)
        ),
    );

    debug.set_match(
        target_idx,
        profile_idx,
        cloud_matrix.get_match(cloud_matrix_row_idx, profile_idx),
    );
    debug.set_insert(
        target_idx,
        profile_idx,
        cloud_matrix.get_insert(cloud_matrix_row_idx, profile_idx),
    );
    debug.set_delete(
        target_idx,
        profile_idx,
        cloud_matrix.get_delete(cloud_matrix_row_idx, profile_idx),
    );
}

#[inline]
pub fn prune_and_scrub(
    bound: &mut CloudBound,
    cloud_matrix: &mut CloudMatrixLinear,
    row_idx: usize,
    alpha: f32,
    beta: f32,
    overall_max: &mut f32,
) {
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
}

pub fn cloud_search_forward(
    profile: &Profile,
    target: &Sequence,
    cloud_matrix: &mut CloudMatrixLinear,
    params: &CloudSearchParams,
    bounds: &mut CloudBoundGroup,
) -> Result<()> {
    let mut debug = DpMatrix::new(target.length, profile.length);

    debug.set_match(params.profile_start, params.target_start, 0.0);
    debug.set_insert(params.profile_start, params.target_start, 0.0);
    debug.set_delete(params.profile_start, params.target_start, 0.0);

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
    cloud_matrix.set_match(first_cloud_matrix_row_idx, params.target_start, 0.0);
    cloud_matrix.set_insert(first_cloud_matrix_row_idx, params.target_start, 0.0);
    cloud_matrix.set_delete(first_cloud_matrix_row_idx, params.target_start, 0.0);

    // the first bound is just the starting cell
    bounds.set(
        first_anti_diagonal_idx,
        params.profile_start,
        params.target_start,
        params.profile_start,
        params.target_start,
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
                &mut debug,
            );
        }
    }
    // main recursion:
    //  - compute the next full anti-diagonal
    //  - prune with max-in-anti-diagonal rule: discard any value < max_in_current_diagonal - alpha (default: alpha = 12)
    //  - prune with x-drop rule: discard anything value < max_in_all_diagonals - beta (default: beta = 20)
    for anti_diagonal_idx in gamma_anti_diagonal_idx..=max_anti_diagonal_idx {
        let previous_bound = bounds.get(anti_diagonal_idx - 1);

        if previous_bound.was_pruned() {
            bounds.max_anti_diagonal_idx = anti_diagonal_idx - 2;
            break;
        }

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

        for (profile_idx, target_idx) in current_bound.anti_diagonal_cell_zip() {
            compute_forward_cell(
                target,
                profile,
                cloud_matrix,
                cloud_matrix_row_idx,
                target_idx,
                profile_idx,
                &mut debug,
            );
        }

        prune_and_scrub(
            current_bound,
            cloud_matrix,
            cloud_matrix_row_idx,
            params.alpha,
            params.beta,
            &mut overall_max_score,
        );
    }

    // let mut debug_out = BufWriter::new(File::create("./cloud-forward.mtx")?);
    // debug.dump(&mut debug_out)?;

    Ok(())
}

pub fn cloud_search_backward(
    profile: &Profile,
    target: &Sequence,
    cloud_matrix: &mut CloudMatrixLinear,
    params: &CloudSearchParams,
    bounds: &mut CloudBoundGroup,
) -> Result<()> {
    let mut debug = DpMatrix::new(target.length, profile.length);

    debug.set_match(params.profile_end, params.target_end, 0.0);
    debug.set_insert(params.profile_end, params.target_end, 0.0);
    debug.set_delete(params.profile_end, params.target_end, 0.0);

    // the highest score we've seen overall
    let mut overall_max_score = -f32::INFINITY;

    let first_anti_diagonal_idx = params.target_end + params.profile_end - 2;
    let gamma_anti_diagonal_idx = first_anti_diagonal_idx - params.gamma;
    let min_anti_diagonal_idx = 0usize;

    let first_cloud_matrix_row_idx = first_anti_diagonal_idx % 3;

    // setting the scores to 0 is like setting
    // the log odds ratio to 1 since log(0) = 1
    cloud_matrix.set_match(first_cloud_matrix_row_idx, params.target_end, 0.0);
    cloud_matrix.set_insert(first_cloud_matrix_row_idx, params.target_end, 0.0);
    cloud_matrix.set_delete(first_cloud_matrix_row_idx, params.target_end, 0.0);

    // the first bound is the ending position
    bounds.set(
        first_anti_diagonal_idx,
        params.target_end.min(target.length - 1),
        params.target_end.min(profile.length - 1),
        params.target_end.min(target.length - 1),
        params.target_end.min(profile.length - 1),
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
                &mut debug,
            );
        }
    }

    // main recursion:
    //  - compute the next full anti-diagonal
    //  - prune with max-in-anti-diagonal rule: discard any value < max_in_current_diagonal - alpha (default: alpha = 12)
    //  - prune with x-drop rule: discard anything value < max_in_all_diagonals - beta (default: beta = 20)
    for anti_diagonal_idx in (min_anti_diagonal_idx..gamma_anti_diagonal_idx).rev() {
        let cloud_matrix_row_idx = anti_diagonal_idx % 3;
        let previous_bound = bounds.get(anti_diagonal_idx + 1);

        if previous_bound.was_pruned() {
            bounds.min_anti_diagonal_idx = anti_diagonal_idx + 2;
            break;
        }

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

        for i in 0..cloud_matrix.data[0].match_vector.len() {
            cloud_matrix.data[cloud_matrix_row_idx].match_vector[i] = -f32::INFINITY;
            cloud_matrix.data[cloud_matrix_row_idx].insert_vector[i] = -f32::INFINITY;
            cloud_matrix.data[cloud_matrix_row_idx].delete_vector[i] = -f32::INFINITY;
        }

        for (target_idx, profile_idx) in current_bound.anti_diagonal_cell_zip() {
            compute_backward_cell(
                target,
                profile,
                cloud_matrix,
                cloud_matrix_row_idx,
                target_idx,
                profile_idx,
                &mut debug,
            );
        }

        prune_and_scrub(
            current_bound,
            cloud_matrix,
            cloud_matrix_row_idx,
            params.alpha,
            params.beta,
            &mut overall_max_score,
        );
    }

    // let mut debug_out = BufWriter::new(File::create("./cloud-backward.mtx")?);
    // debug.dump(&mut debug_out)?;

    Ok(())
}

pub fn cloud_search(
    profile: &Profile,
    target: &Sequence,
    cloud_matrix: &mut CloudMatrixLinear,
    params: &CloudSearchParams,
) -> Result<CloudBoundGroup> {
    let mut forward_bounds = CloudBoundGroup::new(target.length, profile.length);
    cloud_search_forward(profile, target, cloud_matrix, params, &mut forward_bounds)?;

    // let mut forward_json_out = BufWriter::new(File::create("./fwd.json")?);
    // forward_bounds.soda_json(&mut forward_json_out)?;

    cloud_matrix.reuse();
    let mut backward_bounds = CloudBoundGroup::new(target.length, profile.length);
    cloud_search_backward(profile, target, cloud_matrix, params, &mut backward_bounds)?;

    // let mut backward_json_out = BufWriter::new(File::create("./bwd.json")?);
    // backward_bounds.soda_json(&mut backward_json_out)?;

    join_bounds(&mut forward_bounds, &backward_bounds)?;
    forward_bounds.trim_wings();

    // let mut bounds_json_out = BufWriter::new(File::create("./joined.json")?);
    // forward_bounds.soda_json(&mut bounds_json_out)?;

    Ok(forward_bounds)
}
