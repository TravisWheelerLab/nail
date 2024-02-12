use crate::align::structs::{AntiDiagonal, CloudMatrixLinear};
use crate::max_f32;

pub enum PruneStatus {
    FullyPruned,
    PartiallyPruned,
}

pub struct CloudSearchScores {
    pub max_score: f32,
    pub max_score_within: f32,
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
