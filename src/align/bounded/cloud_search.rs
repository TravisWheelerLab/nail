use crate::align::bounded::structs::{CloudBound, CloudMatrix, CloudSearchParams};
use crate::structs::profile::constants::{
    PROFILE_DELETE_TO_DELETE, PROFILE_DELETE_TO_MATCH, PROFILE_INSERT_TO_INSERT,
    PROFILE_INSERT_TO_MATCH, PROFILE_MATCH_TO_DELETE, PROFILE_MATCH_TO_INSERT,
    PROFILE_MATCH_TO_MATCH,
};
use crate::structs::{DpMatrix, Profile, Sequence};
use crate::util::{log_add, PrintMe};
use crate::{log_sum, max_f32};
use anyhow::Result;
use std::fs::File;
use std::io::BufWriter;

#[inline]
pub fn compute_forward_cell(
    profile: &Profile,
    target: &Sequence,
    cloud_matrix: &mut CloudMatrix,
    row_idx: usize,
    target_idx: usize,
    profile_idx: usize,
    debug: &mut DpMatrix,
) {
    let current_target_character = target.digital_bytes[target_idx];

    // match state
    // note: begin to match excluded here
    //
    // in the classic orientation, the match state depends on the
    // cell to the upper left i.e. (row_idx - 1, profile_idx - 1)
    //
    // here, it's the upper left 2-back cell,
    let mut source_row_idx = (row_idx + 1) % 3;
    let mut source_profile_idx = profile_idx - 1;
    cloud_matrix.set_match(
        row_idx,
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
    // in the classic orientation, the insert state depends
    // on the cell above i.e. (row_idx - 1, profile_idx)
    //
    // here, it's actually the same
    source_row_idx = (row_idx + 2) % 3;
    cloud_matrix.set_insert(
        row_idx,
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
    // in the classic orientation, the insert state depends
    // on the cell to the left i.e. (row_idx, profile_idx - 1)
    //
    // here, it's going to be the upper left
    source_row_idx = (row_idx + 2) % 3;
    source_profile_idx = profile_idx - 1;
    cloud_matrix.set_delete(
        row_idx,
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
        cloud_matrix.get_match(row_idx, profile_idx),
    );
    debug.set_insert(
        target_idx,
        profile_idx,
        cloud_matrix.get_insert(row_idx, profile_idx),
    );
    debug.set_delete(
        target_idx,
        profile_idx,
        cloud_matrix.get_delete(row_idx, profile_idx),
    );
}

#[inline]
pub fn compute_backward_cell(
    profile: &Profile,
    target: &Sequence,
    cloud_matrix: &mut CloudMatrix,
    row_idx: usize,
    target_idx: usize,
    profile_idx: usize,
    debug: &mut DpMatrix,
) {
    let current_target_character = target.digital_bytes[target_idx + 1] as usize;

    // match state
    let mut source_row_idx = (row_idx + 1) % 3;
    cloud_matrix.set_match(
        row_idx,
        profile_idx,
        log_sum!(
            cloud_matrix.get_match(source_row_idx, profile_idx + 1)
                + profile.transition_score(PROFILE_MATCH_TO_MATCH, profile_idx)
                + profile.match_score(current_target_character, profile_idx + 1),
            cloud_matrix.get_insert(source_row_idx, profile_idx)
                + profile.transition_score(PROFILE_MATCH_TO_INSERT, profile_idx)
                + profile.insert_score(current_target_character, profile_idx),
            cloud_matrix.get_delete(row_idx, profile_idx + 1)
                + profile.transition_score(PROFILE_MATCH_TO_DELETE, profile_idx)
        ),
    );

    // insert state
    source_row_idx = (row_idx + 2) % 3;
    cloud_matrix.set_insert(
        row_idx,
        profile_idx,
        log_sum!(
            cloud_matrix.get_match(source_row_idx, profile_idx + 1)
                + profile.transition_score(PROFILE_INSERT_TO_MATCH, profile_idx)
                + profile.match_score(current_target_character, profile_idx + 1),
            cloud_matrix.get_insert(source_row_idx, profile_idx)
                + profile.transition_score(PROFILE_INSERT_TO_INSERT, profile_idx)
                + profile.insert_score(current_target_character, profile_idx)
        ),
    );

    // delete state
    source_row_idx = (row_idx + 2) % 3;
    cloud_matrix.set_delete(
        row_idx,
        profile_idx,
        log_sum!(
            cloud_matrix.get_match(source_row_idx, profile_idx + 1)
                + profile.transition_score(PROFILE_DELETE_TO_MATCH, profile_idx)
                + profile.match_score(current_target_character, profile_idx + 1),
            cloud_matrix.get_delete(source_row_idx, profile_idx + 1)
                + profile.transition_score(PROFILE_DELETE_TO_DELETE, profile_idx)
        ),
    );

    debug.set_match(
        target_idx,
        profile_idx,
        cloud_matrix.get_match(row_idx, profile_idx),
    );
    debug.set_insert(
        target_idx,
        profile_idx,
        cloud_matrix.get_insert(row_idx, profile_idx),
    );
    debug.set_delete(
        target_idx,
        profile_idx,
        cloud_matrix.get_delete(row_idx, profile_idx),
    );
}

#[inline]
pub fn prune_and_scrub(
    bound: &mut CloudBound,
    cloud_matrix: &mut CloudMatrix,
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
            bound.left_target_idx -= 1;
            bound.left_profile_idx += 1;
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
            bound.right_target_idx += 1;
            bound.right_profile_idx -= 1;
            cloud_matrix.set_match(row_idx, profile_idx, -f32::INFINITY);
            cloud_matrix.set_insert(row_idx, profile_idx, -f32::INFINITY);
            cloud_matrix.set_delete(row_idx, profile_idx, -f32::INFINITY);
        }
    }
}

pub fn cloud_search_forward(
    profile: &Profile,
    target: &Sequence,
    cloud_matrix: &mut CloudMatrix,
    params: &CloudSearchParams,
) -> Result<Vec<CloudBound>> {
    let mut debug = DpMatrix::new(profile.length, target.length);

    let mut bounds = vec![];

    // setting the scores to 0 is like setting
    // the log odds ratio to 1 since log(0) = 1
    cloud_matrix.set_match(0, params.target_start, 0.0);
    cloud_matrix.set_insert(0, params.target_start, 0.0);
    cloud_matrix.set_delete(0, params.target_start, 0.0);

    debug.set_match(params.profile_start, params.target_start, 0.0);
    debug.set_insert(params.profile_start, params.target_start, 0.0);
    debug.set_delete(params.profile_start, params.target_start, 0.0);

    // initial net: build up to length gamma (default: gamma = 5) anti-diagonal
    //        first anti-diagonal
    //       /   second anti-diagonal
    //      /  /   and so on
    //     /  /  /
    // 0  1  2  3  4  -
    // 1  2  3  4  -  -
    // 2  3  4  -  -  -
    // 3  4  -  -  -  -
    // 4  -  -  -  -  -
    // -  -  -  -  -  -

    // the first bound is the starting position
    bounds.push(CloudBound {
        left_target_idx: params.target_start,
        left_profile_idx: params.profile_start,
        right_target_idx: params.target_start,
        right_profile_idx: params.profile_start,
    });

    // the index of the row in the cloud matrix
    let mut row_idx;
    // the highest score we've seen overall
    let mut overall_max = -f32::INFINITY;

    for anti_diagonal_idx in 1..params.gamma {
        let bound = CloudBound {
            left_target_idx: params.target_start + anti_diagonal_idx,
            left_profile_idx: params.profile_start,
            right_target_idx: params.target_start,
            right_profile_idx: params.profile_start + anti_diagonal_idx,
        };

        row_idx = anti_diagonal_idx % 3;

        for (target_idx, profile_idx) in bound.anti_diagonal() {
            compute_forward_cell(
                profile,
                target,
                cloud_matrix,
                row_idx,
                target_idx,
                profile_idx,
                &mut debug,
            );
        }
        bounds.push(bound);
    }

    // main recursion:
    //  - compute the next full anti-diagonal
    //  - prune with max-in-anti-diagonal rule: discard any value < max_in_current_diagonal - alpha (default: alpha = 12)
    //  - prune with x-drop rule: discard anything value < max_in_all_diagonals - beta (default: beta = 20)
    // TODO: these loop bounds are wrong for anything other than starting at (1, 1)
    for anti_diagonal_idx in params.gamma..(target.length + profile.length) {
        row_idx = anti_diagonal_idx % 3;

        let previous_bound = &bounds[anti_diagonal_idx - 1];

        // TODO: there's a slight chance that looping a no-op might
        //       be faster than keeping this if statement here
        if previous_bound.was_pruned() {
            break;
        }

        // TODO: maybe it would be faster to allocate all
        //      of the bounds at start, then adjust them?
        let mut bound = previous_bound.clone();

        // edge guards
        // TODO: I feel like this could be vastly improved,
        //       maybe without any conditionals
        if bound.left_target_idx < target.length {
            // if we haven't hit the end of the target (bottom row),
            // we'll plan to move the left bound down 1
            bound.left_target_idx += 1;
        } else {
            // otherwise we'll move to the right 1
            bound.left_profile_idx += 1;
        }

        if bound.right_profile_idx < profile.length {
            // if we haven't hit the end of the profile (last column),
            // we'll plan to move the right bound to the right 1
            bound.right_profile_idx += 1;
        } else {
            // otherwise we'll move down 1
            bound.right_target_idx += 1;
        }

        for (target_idx, profile_idx) in bound.anti_diagonal() {
            compute_forward_cell(
                profile,
                target,
                cloud_matrix,
                row_idx,
                target_idx,
                profile_idx,
                &mut debug,
            );
        }

        prune_and_scrub(
            &mut bound,
            cloud_matrix,
            row_idx,
            params.alpha,
            params.beta,
            &mut overall_max,
        );

        bounds.push(bound);
    }

    let mut debug_out = BufWriter::new(File::create("./cloud-forward.mtx")?);
    debug.dump(&mut debug_out)?;

    Ok(bounds)
}

pub fn cloud_search_backward(
    profile: &Profile,
    target: &Sequence,
    cloud_matrix: &mut CloudMatrix,
    params: &CloudSearchParams,
) -> Result<Vec<CloudBound>> {
    let mut debug = DpMatrix::new(profile.length, target.length);

    let mut bounds = vec![];

    // setting the scores to 0 is like setting
    // the log odds ratio to 1 since log(0) = 1
    cloud_matrix.set_match(0, params.target_end, 0.0);
    cloud_matrix.set_insert(0, params.target_end, 0.0);
    cloud_matrix.set_delete(0, params.target_end, 0.0);

    debug.set_match(params.profile_end, params.target_end, 0.0);
    debug.set_insert(params.profile_end, params.target_end, 0.0);
    debug.set_delete(params.profile_end, params.target_end, 0.0);

    // initial net: build up to length gamma (default: gamma = 5) anti-diagonal
    // -  -  -  -  -  -
    // -  -  -  -  -  4
    // -  -  -  -  4  3
    // -  -  -  4  3  2
    // -  -  4  3  2  1
    // -  4  3  2  1  0

    // the first bound is the ending position
    bounds.push(CloudBound {
        left_target_idx: params.target_end,
        left_profile_idx: params.profile_end,
        right_target_idx: params.target_end,
        right_profile_idx: params.profile_end,
    });

    // the index of the row in the cloud matrix
    let mut row_idx;
    // the highest score we've seen overall
    let mut overall_max = -f32::INFINITY;

    for anti_diagonal_idx in 1..params.gamma {
        let bound = CloudBound {
            left_target_idx: params.target_end,
            left_profile_idx: params.profile_end - anti_diagonal_idx,
            right_target_idx: params.target_end - anti_diagonal_idx,
            right_profile_idx: params.profile_end,
        };
        row_idx = anti_diagonal_idx % 3;

        for (target_idx, profile_idx) in bound.anti_diagonal() {
            compute_backward_cell(
                profile,
                target,
                cloud_matrix,
                row_idx,
                target_idx,
                profile_idx,
                &mut debug,
            );
        }
        bounds.push(bound);
    }

    // main recursion:
    //  - compute the next full anti-diagonal
    //  - prune with max-in-anti-diagonal rule: discard any value < max_in_current_diagonal - alpha (default: alpha = 12)
    //  - prune with x-drop rule: discard anything value < max_in_all_diagonals - beta (default: beta = 20)
    // TODO: these loop bounds are wrong for anything other than starting at (1, 1)

    for anti_diagonal_idx in params.gamma..(target.length + profile.length) {
        row_idx = anti_diagonal_idx % 3;

        let previous_bound = &bounds[anti_diagonal_idx - 1];

        if previous_bound.was_pruned() {
            break;
        }

        let mut bound = previous_bound.clone();

        // edge guards
        if bound.left_profile_idx > 1 {
            // if we haven't hit the start of the target (top row),
            // we'll plan to move the left bound up
            bound.left_profile_idx -= 1;
        } else {
            // otherwise we'll move to the left
            bound.left_target_idx -= 1;
        }

        if bound.right_target_idx > 1 {
            // if we haven't hit the start of the profile (first column),
            // we'll plan to move the right bound to the left
            bound.right_target_idx -= 1;
        } else {
            // otherwise we'll move up
            bound.right_profile_idx -= 1;
        }

        for i in 0..cloud_matrix.data[0].match_vector.len() {
            cloud_matrix.data[row_idx].match_vector[i] = -f32::INFINITY;
            cloud_matrix.data[row_idx].insert_vector[i] = -f32::INFINITY;
            cloud_matrix.data[row_idx].delete_vector[i] = -f32::INFINITY;
        }

        for (target_idx, profile_idx) in bound.anti_diagonal() {
            compute_backward_cell(
                profile,
                target,
                cloud_matrix,
                row_idx,
                target_idx,
                profile_idx,
                &mut debug,
            );
        }

        prune_and_scrub(
            &mut bound,
            cloud_matrix,
            row_idx,
            params.alpha,
            params.beta,
            &mut overall_max,
        );

        bounds.push(bound);
    }

    let mut debug_out = BufWriter::new(File::create("./cloud-backward.mtx")?);
    debug.dump(&mut debug_out)?;

    Ok(bounds)
}

pub fn cloud_search(
    profile: &Profile,
    target: &Sequence,
    cloud_matrix: &mut CloudMatrix,
    params: &CloudSearchParams,
) -> Result<Vec<CloudBound>> {
    let forward_bounds = cloud_search_forward(profile, target, cloud_matrix, params);
    cloud_matrix.reuse();
    let backward_bounds = cloud_search_backward(profile, target, cloud_matrix, params);

    let bounds = vec![];

    Ok(bounds)
}
