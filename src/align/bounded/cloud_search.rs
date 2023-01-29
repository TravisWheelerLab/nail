use crate::align::bounded::structs::{Bound, CloudMatrix};
use crate::log_sum;
use crate::structs::profile::constants::{
    PROFILE_DELETE_TO_DELETE, PROFILE_DELETE_TO_MATCH, PROFILE_INSERT_TO_INSERT,
    PROFILE_INSERT_TO_MATCH, PROFILE_MATCH_TO_DELETE, PROFILE_MATCH_TO_INSERT,
    PROFILE_MATCH_TO_MATCH,
};
use crate::structs::{DpMatrix, Profile, Sequence};
use crate::util::{log_sum, PrintMe};
use anyhow::Result;
use std::fs::File;
use std::io::BufWriter;
use std::iter::zip;

#[inline]
pub fn compute_forward_cell(
    profile: &Profile,
    target: &Sequence,
    cloud_matrix: &mut CloudMatrix,
    target_idx: usize,
    profile_idx: usize,
    anti_diagonal_idx: usize,
    cell_idx: usize,
    debug: &mut DpMatrix,
) {
    let current_target_character = target.digital_bytes[target_idx];

    let row_idx = anti_diagonal_idx % 3;
    // offset by one for the left buffer
    let col_idx = cell_idx + 1;

    // match state
    // note: begin to match excluded here
    //
    // in the classic orientation, the match state depends on the
    // cell to the upper left i.e. (row_idx - 1, col_idx - 1)
    //
    // here, it's the upper left 2-back cell,
    // i.e. ((row_idx - 2) % 3, col_idx - 1),
    let mut source_row_idx = (row_idx + 1) % 3;
    let mut source_col_idx = (col_idx + 2) % 3;
    cloud_matrix.set_match(
        row_idx,
        col_idx,
        log_sum!(
            cloud_matrix.get_match(source_row_idx, source_col_idx)
                + profile.transition_score(PROFILE_MATCH_TO_MATCH, profile_idx - 1),
            cloud_matrix.get_insert(source_row_idx, source_col_idx)
                + profile.transition_score(PROFILE_INSERT_TO_MATCH, profile_idx - 1),
            cloud_matrix.get_delete(source_row_idx, source_col_idx)
                + profile.transition_score(PROFILE_DELETE_TO_MATCH, profile_idx - 1)
        ) + profile.match_score(current_target_character as usize, profile_idx),
    );

    // insert state
    //
    // in the classic orientation, the insert state depends
    // on the cell above i.e. (row_idx - 1, col_idx)
    //
    // here, it's actually the same
    source_row_idx = (row_idx + 2) % 3;
    source_col_idx = col_idx;
    cloud_matrix.set_insert(
        row_idx,
        col_idx,
        log_sum(
            cloud_matrix.get_match(source_row_idx, source_col_idx)
                + profile.transition_score(PROFILE_MATCH_TO_INSERT, profile_idx),
            cloud_matrix.get_insert(source_row_idx, source_col_idx)
                + profile.transition_score(PROFILE_INSERT_TO_INSERT, profile_idx),
        ) + profile.insert_score(current_target_character as usize, profile_idx),
    );

    // delete state
    //
    // in the classic orientation, the insert state depends
    // on the cell to the left i.e. (row_idx, col_idx - 1)
    //
    // here, it's going to be the upper left
    // i.e. (row_idx - 1, col_idx - 1)
    source_row_idx = (row_idx + 2) % 3;
    source_col_idx = (col_idx + 2) % 3;
    cloud_matrix.set_delete(
        row_idx,
        col_idx,
        log_sum(
            cloud_matrix.get_match(source_row_idx, source_col_idx)
                + profile.transition_score(PROFILE_MATCH_TO_DELETE, profile_idx - 1),
            cloud_matrix.get_delete(source_row_idx, source_col_idx)
                + profile.transition_score(PROFILE_DELETE_TO_DELETE, profile_idx - 1),
        ),
    );

    debug.set_match(
        target_idx,
        profile_idx,
        cloud_matrix.get_match(row_idx, col_idx),
    );
    debug.set_insert(
        target_idx,
        profile_idx,
        cloud_matrix.get_insert(row_idx, col_idx),
    );
    debug.set_delete(
        target_idx,
        profile_idx,
        cloud_matrix.get_delete(row_idx, col_idx),
    );
}

pub fn cloud_search(
    profile: &Profile,
    target: &Sequence,
    cloud_matrix: &mut CloudMatrix,
) -> Result<Vec<Bound>> {
    let mut debug = DpMatrix::new(profile.length, target.length);

    let bounds = vec![];

    // TODO: parameterize these
    let target_start: usize = 1;
    let profile_start: usize = 1;
    let gamma: usize = 5;
    let alpha: f32 = 12.0;
    let beta: f32 = 20.0;

    // instead of initializing a zero row, we just want to set the starting cell
    // start at col_idx = 1 because col_idx = 0 is the buffer of -inf
    cloud_matrix.set_match(0, 1, 0.0);
    cloud_matrix.set_insert(0, 1, 0.0);
    cloud_matrix.set_delete(0, 1, 0.0);

    debug.set_match(1, 1, 0.0);
    debug.set_insert(1, 1, 0.0);
    debug.set_delete(1, 1, 0.0);

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

    let mut left_bound = Bound {
        target_idx: target_start,
        profile_idx: profile_start,
    };

    let mut right_bound = Bound {
        target_idx: target_start,
        profile_idx: profile_start,
    };

    // the index of the current anti-diagonal
    let mut anti_diagonal_idx = 0;

    // let mut overall_max = -f32::INFINITY;

    for _ in 1..gamma {
        anti_diagonal_idx += 1;
        left_bound.target_idx += 1;
        right_bound.profile_idx += 1;

        let target_range = (right_bound.target_idx..=left_bound.target_idx).rev();
        let profile_range = left_bound.profile_idx..=right_bound.profile_idx;

        for (cell_idx, (target_idx, profile_idx)) in zip(target_range, profile_range).enumerate() {
            compute_forward_cell(
                profile,
                target,
                cloud_matrix,
                target_idx,
                profile_idx,
                anti_diagonal_idx,
                cell_idx,
                &mut debug,
            );
        }
    }

    // main recursion:
    // - compute the next full anti-diagonal (should be previous anti-diagonal length + 1)
    // - prune with max-in-anti-diagonal rule: discard any value < max_in_current_diagonal - alpha (default: alpha = 12)
    // - prune with x-drop rule: discard anything value < max_in_all_diagonals - beta (default: beta = 20)
    loop {
        if left_bound.target_idx < right_bound.target_idx {
            break;
        }

        if left_bound.target_idx == right_bound.target_idx
            && left_bound.target_idx == target.length
            && left_bound.profile_idx == right_bound.profile_idx
            && left_bound.profile_idx == profile.length
        {
            break;
        }

        let mut current_max = -f32::INFINITY;

        // edge guards
        if left_bound.target_idx < target.length {
            left_bound.target_idx += 1;
        } else {
            left_bound.profile_idx += 1;
        }

        if right_bound.profile_idx < profile.length {
            right_bound.profile_idx += 1;
        } else {
            right_bound.target_idx += 1;
        }

        let num_cells = right_bound.profile_idx - left_bound.profile_idx + 1;
        let target_range = (right_bound.target_idx..=left_bound.target_idx).rev();
        let profile_range = left_bound.profile_idx..=right_bound.profile_idx;

        for (cell_idx, (target_idx, profile_idx)) in zip(target_range, profile_range).enumerate() {
            compute_forward_cell(
                profile,
                target,
                cloud_matrix,
                target_idx,
                profile_idx,
                anti_diagonal_idx,
                cell_idx,
                &mut debug,
            );
        }

        let match_values = &cloud_matrix.data[anti_diagonal_idx % 3].match_vector[1..=num_cells];

        for value in match_values.iter() {
            current_max = current_max.max(*value);
        }

        'inner: for value in match_values.iter() {
            if *value <= current_max - alpha {
                // println!("L: {} pruned", value);
                left_bound.target_idx -= 1;
                left_bound.profile_idx += 1;
            } else {
                break 'inner;
            }
        }

        'inner: for value in match_values.iter().rev() {
            if *value <= current_max - alpha {
                // println!("R: {} pruned", value);
                right_bound.target_idx += 1;
                right_bound.profile_idx -= 1;
            } else {
                break 'inner;
            }
        }
    }

    let mut debug_out = BufWriter::new(File::create("./cloud-debug.mtx")?);
    debug.dump(&mut debug_out)?;
    Ok(bounds)
}
