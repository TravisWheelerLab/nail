use crate::align::structs::{DpMatrix, RowBounds};
use crate::log_sum;
use crate::structs::{Profile, Sequence};
use crate::util::log_add;

pub fn backward(
    profile: &Profile,
    target: &Sequence,
    dp_matrix: &mut impl DpMatrix,
    row_bounds: &RowBounds,
) {
    dp_matrix.set_special(row_bounds.seq_end, Profile::B_IDX, -f32::INFINITY);
    dp_matrix.set_special(row_bounds.seq_end, Profile::N_IDX, -f32::INFINITY);

    dp_matrix.set_special(
        row_bounds.seq_end,
        Profile::C_IDX,
        profile.special_transition_score(Profile::C_IDX, Profile::SPECIAL_MOVE_IDX),
    );

    dp_matrix.set_special(
        row_bounds.seq_end,
        Profile::E_IDX,
        profile.special_transition_score(Profile::C_IDX, Profile::SPECIAL_MOVE_IDX)
            + profile.special_transition_score(Profile::E_IDX, Profile::SPECIAL_MOVE_IDX),
    );

    let profile_start_on_last_row = row_bounds.left_row_bounds[row_bounds.seq_end];
    let profile_end_on_last_row = row_bounds.right_row_bounds[row_bounds.seq_end];

    // C: dp matrix cell
    // L: cell in the last row
    // *: the last cell on the last row
    //
    //     bounded                       full
    //  C  C  C  C  -  -            F  F  F  F  F  F
    //  C  C  C  C  -  -            F  F  F  F  F  F
    //  C  C  C  C  -  -   <-vs->   F  F  F  F  F  F
    //  L  L  L  *  -  -            F  F  F  F  F  F
    //  -  - -  -  -  -             L  L  L  L  L  *

    dp_matrix.set_match(
        row_bounds.seq_end,
        profile_end_on_last_row,
        dp_matrix.get_special(row_bounds.seq_end, Profile::E_IDX),
    );

    dp_matrix.set_delete(
        row_bounds.seq_end,
        profile_end_on_last_row,
        dp_matrix.get_special(row_bounds.seq_end, Profile::E_IDX),
    );

    dp_matrix.set_insert(row_bounds.seq_end, profile_end_on_last_row, -f32::INFINITY);

    // this loops over the last row, setting all of the <L> cells, excluding the last cell <*>
    for profile_idx in (profile_start_on_last_row..profile_end_on_last_row).rev() {
        dp_matrix.set_match(
            row_bounds.seq_end,
            profile_idx,
            log_sum!(
                dp_matrix.get_special(row_bounds.seq_end, Profile::E_IDX),
                dp_matrix.get_delete(row_bounds.seq_end, profile_idx + 1)
                    + profile.transition_score(Profile::M_D_IDX, profile_idx)
            ),
        );

        dp_matrix.set_insert(row_bounds.seq_end, profile_idx, -f32::INFINITY);
        dp_matrix.set_delete(
            row_bounds.seq_end,
            profile_idx,
            log_sum!(
                dp_matrix.get_special(row_bounds.seq_end, Profile::E_IDX),
                dp_matrix.get_delete(row_bounds.seq_end, profile_idx + 1)
                    + profile.transition_score(Profile::D_D_IDX, profile_idx)
            ),
        );
    }

    // main recursion
    for target_idx in (row_bounds.seq_start..row_bounds.seq_end).rev() {
        let current_residue = target.digital_bytes[target_idx + 1] as usize;
        let profile_start_on_current_row = row_bounds.left_row_bounds[target_idx];
        let profile_end_on_current_row = row_bounds.right_row_bounds[target_idx];

        // TODO: I don't think we need to do this as long as the B state is initialized with -inf
        //       if we do need to do this for some reason, we need to change the bounds of the next loop
        dp_matrix.set_special(
            target_idx,
            Profile::B_IDX,
            dp_matrix.get_match(target_idx + 1, profile_start_on_current_row)
                + profile.transition_score(Profile::B_M_IDX, profile_start_on_current_row - 1)
                + profile.match_score(current_residue, profile_start_on_current_row),
        );

        // this loops over the cells in the current row, incrementally adding to the B state
        // it depends on the match scores in the next row being computed first
        for profile_idx in (profile_start_on_current_row + 1)..=profile_end_on_current_row {
            dp_matrix.set_special(
                target_idx,
                Profile::B_IDX,
                log_sum!(
                    dp_matrix.get_special(target_idx, Profile::B_IDX),
                    dp_matrix.get_match(target_idx + 1, profile_idx)
                        + profile.transition_score(Profile::B_M_IDX, profile_idx - 1)
                        + profile.match_score(current_residue, profile_idx)
                ),
            );
        }

        dp_matrix.set_special(
            target_idx,
            Profile::C_IDX,
            dp_matrix.get_special(target_idx + 1, Profile::C_IDX)
                + profile.special_transition_score(Profile::C_IDX, Profile::SPECIAL_LOOP_IDX),
        );

        dp_matrix.set_special(
            target_idx,
            Profile::E_IDX,
            log_sum!(
                dp_matrix.get_special(target_idx, Profile::C_IDX)
                    + profile.special_transition_score(Profile::E_IDX, Profile::SPECIAL_MOVE_IDX)
            ),
        );

        dp_matrix.set_special(
            target_idx,
            Profile::N_IDX,
            log_sum!(
                dp_matrix.get_special(target_idx + 1, Profile::N_IDX)
                    + profile.special_transition_score(Profile::N_IDX, Profile::SPECIAL_LOOP_IDX),
                dp_matrix.get_special(target_idx, Profile::B_IDX)
                    + profile.special_transition_score(Profile::N_IDX, Profile::SPECIAL_MOVE_IDX)
            ),
        );

        dp_matrix.set_match(
            target_idx,
            profile_end_on_current_row,
            dp_matrix.get_special(target_idx, Profile::E_IDX),
        );

        dp_matrix.set_insert(target_idx, profile_end_on_current_row, -f32::INFINITY);

        dp_matrix.set_delete(
            target_idx,
            profile_end_on_current_row,
            dp_matrix.get_special(target_idx, Profile::E_IDX),
        );

        for profile_idx in (profile_start_on_current_row..profile_end_on_current_row).rev() {
            dp_matrix.set_match(
                target_idx,
                profile_idx,
                log_sum!(
                    dp_matrix.get_match(target_idx + 1, profile_idx + 1)
                        + profile.transition_score(Profile::M_M_IDX, profile_idx)
                        + profile.match_score(current_residue, profile_idx + 1),
                    dp_matrix.get_insert(target_idx + 1, profile_idx)
                        + profile.transition_score(Profile::M_I_IDX, profile_idx)
                        + profile.insert_score(current_residue, profile_idx),
                    dp_matrix.get_special(target_idx, Profile::E_IDX),
                    dp_matrix.get_delete(target_idx, profile_idx + 1)
                        + profile.transition_score(Profile::M_D_IDX, profile_idx)
                ),
            );

            dp_matrix.set_insert(
                target_idx,
                profile_idx,
                log_sum!(
                    dp_matrix.get_match(target_idx + 1, profile_idx + 1)
                        + profile.transition_score(Profile::I_M_IDX, profile_idx)
                        + profile.match_score(current_residue, profile_idx + 1),
                    dp_matrix.get_insert(target_idx + 1, profile_idx)
                        + profile.transition_score(Profile::I_I_IDX, profile_idx)
                        + profile.insert_score(current_residue, profile_idx)
                ),
            );

            dp_matrix.set_delete(
                target_idx,
                profile_idx,
                log_sum!(
                    dp_matrix.get_match(target_idx + 1, profile_idx + 1)
                        + profile.transition_score(Profile::D_M_IDX, profile_idx)
                        + profile.match_score(current_residue, profile_idx + 1),
                    dp_matrix.get_delete(target_idx, profile_idx + 1)
                        + profile.transition_score(Profile::D_D_IDX, profile_idx),
                    dp_matrix.get_special(target_idx, Profile::E_IDX)
                ),
            );
        }
    }

    let first_target_character = target.digital_bytes[row_bounds.seq_start] as usize;

    let profile_start_in_first_row = row_bounds.left_row_bounds[row_bounds.seq_start];
    let profile_end_in_first_row = row_bounds.right_row_bounds[row_bounds.seq_start];

    dp_matrix.set_special(
        row_bounds.seq_start - 1,
        Profile::B_IDX,
        dp_matrix.get_match(row_bounds.seq_start, profile_start_in_first_row)
            + profile.transition_score(Profile::B_M_IDX, 0)
            + profile.match_score(first_target_character, 1),
    );

    for profile_idx in (profile_start_in_first_row + 1)..=profile_end_in_first_row {
        dp_matrix.set_special(
            row_bounds.seq_start - 1,
            Profile::B_IDX,
            log_sum!(
                dp_matrix.get_special(row_bounds.seq_start - 1, Profile::B_IDX),
                dp_matrix.get_match(row_bounds.seq_start, profile_idx)
                    + profile.transition_score(Profile::B_M_IDX, profile_idx - 1)
                    + profile.match_score(first_target_character, profile_idx)
            ),
        );
    }

    // dp_matrix.set_special(
    //     row_bounds.target_start - 1,
    //     Profile::SPECIAL_J_IDX,
    //     -f32::INFINITY,
    // );
    dp_matrix.set_special(row_bounds.seq_start - 1, Profile::C_IDX, -f32::INFINITY);
    dp_matrix.set_special(row_bounds.seq_start - 1, Profile::E_IDX, -f32::INFINITY);
    dp_matrix.set_special(
        row_bounds.seq_start - 1,
        Profile::N_IDX,
        log_sum!(
            dp_matrix.get_special(row_bounds.seq_start, Profile::N_IDX)
                + profile.special_transition_score(Profile::N_IDX, Profile::SPECIAL_LOOP_IDX),
            dp_matrix.get_special(row_bounds.seq_start - 1, Profile::B_IDX)
                + profile.special_transition_score(Profile::N_IDX, Profile::SPECIAL_MOVE_IDX)
        ),
    );
    for profile_idx in (profile_start_in_first_row..=profile_end_in_first_row).rev() {
        dp_matrix.set_match(row_bounds.seq_start - 1, profile_idx, -f32::INFINITY);
        dp_matrix.set_insert(row_bounds.seq_start - 1, profile_idx, -f32::INFINITY);
        dp_matrix.set_delete(row_bounds.seq_start - 1, profile_idx, -f32::INFINITY);
    }
}
