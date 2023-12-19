use crate::align::structs::{DpMatrix, RowBounds};
use crate::log_sum;
use crate::structs::{Profile, Sequence};
use crate::util::log_add;

pub fn backward_bounded(
    profile: &Profile,
    target: &Sequence,
    dp_matrix: &mut impl DpMatrix,
    row_bounds: &RowBounds,
) {
    let end_score: f32 = 0.0;
    //  M   s_D
    //  s_I s_M
    //
    //  I    -
    //  s_I s_M
    //
    //  D   s_D
    //  -   s_M
    //

    // C:   a cell in the cloud
    // L:   a cell in the last *row* of the cloud
    // [p]: a padding cell
    //
    //  .........................
    //  .........................
    //  ....  C  C  C  C  C  C  C [p]
    //  ....  C  C  C  C  C  C [p][p] <-  we need different amounts of
    //  ....  C  C  C  C  C  C [p]    <-  padding depending on whether
    //  ....  C  C  C  C  C [p][p]    <-  or not the rows are offset
    //  ....  L  L  L  L [p][p]
    //        ^  ^  ^  ^
    //       the last row defined in the row bounds
    //       is initialized with special conditions

    dp_matrix.set_special(
        row_bounds.target_end,
        Profile::SPECIAL_J_IDX,
        -f32::INFINITY,
    );
    dp_matrix.set_special(
        row_bounds.target_end,
        Profile::SPECIAL_B_IDX,
        -f32::INFINITY,
    );
    dp_matrix.set_special(
        row_bounds.target_end,
        Profile::SPECIAL_N_IDX,
        -f32::INFINITY,
    );

    dp_matrix.set_special(
        row_bounds.target_end,
        Profile::SPECIAL_C_IDX,
        profile.special_transition_score(Profile::SPECIAL_C_IDX, Profile::SPECIAL_MOVE_IDX),
    );

    dp_matrix.set_special(
        row_bounds.target_end,
        Profile::SPECIAL_E_IDX,
        profile.special_transition_score(Profile::SPECIAL_C_IDX, Profile::SPECIAL_MOVE_IDX)
            + profile.special_transition_score(Profile::SPECIAL_E_IDX, Profile::SPECIAL_MOVE_IDX),
    );

    let profile_start_on_last_row = row_bounds.left_row_bounds[row_bounds.target_end];
    let profile_end_on_last_row = row_bounds.right_row_bounds[row_bounds.target_end];

    // C: dp matrix cell
    // L: cell in the last row
    // *: the last cell on the last row
    //
    //     bounded                       full
    //  C  C  C  C  C  -            F  F  F  F  F  F
    //  C  C  C  C  -  -            F  F  F  F  F  F
    //  C  C  C  C  -  -   <-vs->   F  F  F  F  F  F
    //  C  C  C  -  -  -            F  F  F  F  F  F
    //  L  *  -  -  -  -            L  L  L  L  L  *

    dp_matrix.set_match(
        row_bounds.target_end,
        profile_end_on_last_row,
        dp_matrix.get_special(row_bounds.target_end, Profile::SPECIAL_E_IDX),
    );

    dp_matrix.set_delete(
        row_bounds.target_end,
        profile_end_on_last_row,
        dp_matrix.get_special(row_bounds.target_end, Profile::SPECIAL_E_IDX),
    );

    dp_matrix.set_insert(
        row_bounds.target_end,
        profile_end_on_last_row,
        -f32::INFINITY,
    );

    // C: dp matrix cell
    // L: cell in the last row
    // *: the last cell on the last row
    //
    //     bounded                       full
    //  C  C  C  C  C  -            F  F  F  F  F  F
    //  C  C  C  C  -  -   <-vs->   F  F  F  F  F  F
    //  C  C  C  -  -  -            F  F  F  F  F  F
    //  L  *  -  -  -  -            F  F  F  F  F  F
    //  -  -  -  -  -  -            L  L  L  L  L  *
    //
    // this loops over the last row, setting all of the <L> cells, excluding the last cell <*>
    for profile_idx in (profile_start_on_last_row..profile_end_on_last_row).rev() {
        dp_matrix.set_match(
            row_bounds.target_end,
            profile_idx,
            log_sum!(
                dp_matrix.get_special(row_bounds.target_end, Profile::SPECIAL_E_IDX) + end_score,
                dp_matrix.get_delete(row_bounds.target_end, profile_idx + 1)
                    + profile.transition_score(Profile::MATCH_TO_DELETE_IDX, profile_idx)
            ),
        );

        dp_matrix.set_insert(row_bounds.target_end, profile_idx, -f32::INFINITY);
        dp_matrix.set_delete(
            row_bounds.target_end,
            profile_idx,
            log_sum!(
                dp_matrix.get_special(row_bounds.target_end, Profile::SPECIAL_E_IDX) + end_score,
                dp_matrix.get_delete(row_bounds.target_end, profile_idx + 1)
                    + profile.transition_score(Profile::DELETE_TO_DELETE_IDX, profile_idx)
            ),
        );
    }

    // main recursion
    for target_idx in (row_bounds.target_start..row_bounds.target_end).rev() {
        let current_residue = target.digital_bytes[target_idx + 1] as usize;
        //            Backward matrix             B state
        // ... .   ...  .  .  .  .  .  .  .  .       .
        // ... .   ...  .  .  .  .  .  .  .  .       .
        // ... .   ...  .  .  .  .  .  .  .  .       .
        // ... C   ...  C  C  C  C  C  C  C  -       -
        // ... C   ...  C  C  C  C  C  C  C  -       -
        // ... C   ...  C  C  C  C  C  C  C  -       -
        // ... C   ...  C  C  C  C  C  C  -  -      B_t = M_p * tsc(B->M) * msc(T_t) // init
        //    M_p  ...  C  C  C  C  C  -  -  -      B_p // the previous B state
        //              .  .  .  .  .  .  .  .       .
        //              .  .  .  .  .  .  .  .       .
        //              .  .  .  .  .  .  .  .       .
        //              L  L  L  L  -  -  -  -      B_L  // the first B state
        //              -  -  -  -  -  -  -  -
        //

        let profile_start_on_current_row = row_bounds.left_row_bounds[target_idx];
        let profile_end_on_current_row = row_bounds.right_row_bounds[target_idx];

        // TODO: I don't think we need to do this as long as the B state is initialized with -inf
        //       if we do need to do this for some reason, we need to change the bounds of the next loop
        dp_matrix.set_special(
            target_idx,
            Profile::SPECIAL_B_IDX,
            dp_matrix.get_match(target_idx + 1, profile_start_on_current_row)
                + profile.transition_score(
                    Profile::BEGIN_TO_MATCH_IDX,
                    profile_start_on_current_row - 1,
                )
                + profile.match_score(current_residue, profile_start_on_current_row),
        );

        // this loops over the cells in the current row, incrementally adding to the B state
        // it depends on the match scores in the next row being computed first
        for profile_idx in (profile_start_on_current_row + 1)..=profile_end_on_current_row {
            dp_matrix.set_special(
                target_idx,
                Profile::SPECIAL_B_IDX,
                log_sum!(
                    dp_matrix.get_special(target_idx, Profile::SPECIAL_B_IDX),
                    dp_matrix.get_match(target_idx + 1, profile_idx)
                        + profile.transition_score(Profile::BEGIN_TO_MATCH_IDX, profile_idx - 1)
                        + profile.match_score(current_residue, profile_idx)
                ),
            );
        }

        dp_matrix.set_special(
            target_idx,
            Profile::SPECIAL_J_IDX,
            log_sum!(
                dp_matrix.get_special(target_idx + 1, Profile::SPECIAL_J_IDX)
                    + profile.special_transition_score(
                        Profile::SPECIAL_J_IDX,
                        Profile::SPECIAL_LOOP_IDX
                    ),
                dp_matrix.get_special(target_idx, Profile::SPECIAL_B_IDX)
                    + profile.special_transition_score(
                        Profile::SPECIAL_J_IDX,
                        Profile::SPECIAL_MOVE_IDX
                    )
            ),
        );

        dp_matrix.set_special(
            target_idx,
            Profile::SPECIAL_C_IDX,
            dp_matrix.get_special(target_idx + 1, Profile::SPECIAL_C_IDX)
                + profile
                    .special_transition_score(Profile::SPECIAL_C_IDX, Profile::SPECIAL_LOOP_IDX),
        );

        dp_matrix.set_special(
            target_idx,
            Profile::SPECIAL_E_IDX,
            log_sum!(
                dp_matrix.get_special(target_idx, Profile::SPECIAL_J_IDX)
                    + profile.special_transition_score(
                        Profile::SPECIAL_E_IDX,
                        Profile::SPECIAL_LOOP_IDX
                    ),
                dp_matrix.get_special(target_idx, Profile::SPECIAL_C_IDX)
                    + profile.special_transition_score(
                        Profile::SPECIAL_E_IDX,
                        Profile::SPECIAL_MOVE_IDX
                    )
            ),
        );

        dp_matrix.set_special(
            target_idx,
            Profile::SPECIAL_N_IDX,
            log_sum!(
                dp_matrix.get_special(target_idx + 1, Profile::SPECIAL_N_IDX)
                    + profile.special_transition_score(
                        Profile::SPECIAL_N_IDX,
                        Profile::SPECIAL_LOOP_IDX
                    ),
                dp_matrix.get_special(target_idx, Profile::SPECIAL_B_IDX)
                    + profile.special_transition_score(
                        Profile::SPECIAL_N_IDX,
                        Profile::SPECIAL_MOVE_IDX
                    )
            ),
        );

        dp_matrix.set_match(
            target_idx,
            profile_end_on_current_row,
            dp_matrix.get_special(target_idx, Profile::SPECIAL_E_IDX),
        );

        dp_matrix.set_insert(target_idx, profile_end_on_current_row, -f32::INFINITY);

        dp_matrix.set_delete(
            target_idx,
            profile_end_on_current_row,
            dp_matrix.get_special(target_idx, Profile::SPECIAL_E_IDX),
        );

        for profile_idx in (profile_start_on_current_row..profile_end_on_current_row).rev() {
            dp_matrix.set_match(
                target_idx,
                profile_idx,
                log_sum!(
                    dp_matrix.get_match(target_idx + 1, profile_idx + 1)
                        + profile.transition_score(Profile::MATCH_TO_MATCH_IDX, profile_idx)
                        + profile.match_score(current_residue, profile_idx + 1),
                    dp_matrix.get_insert(target_idx + 1, profile_idx)
                        + profile.transition_score(Profile::MATCH_TO_INSERT_IDX, profile_idx)
                        + profile.insert_score(current_residue, profile_idx),
                    dp_matrix.get_special(target_idx, Profile::SPECIAL_E_IDX) + end_score,
                    dp_matrix.get_delete(target_idx, profile_idx + 1)
                        + profile.transition_score(Profile::MATCH_TO_DELETE_IDX, profile_idx)
                ),
            );

            dp_matrix.set_insert(
                target_idx,
                profile_idx,
                log_sum!(
                    dp_matrix.get_match(target_idx + 1, profile_idx + 1)
                        + profile.transition_score(Profile::INSERT_TO_MATCH_IDX, profile_idx)
                        + profile.match_score(current_residue, profile_idx + 1),
                    dp_matrix.get_insert(target_idx + 1, profile_idx)
                        + profile.transition_score(Profile::INSERT_TO_INSERT_IDX, profile_idx)
                        + profile.insert_score(current_residue, profile_idx)
                ),
            );

            dp_matrix.set_delete(
                target_idx,
                profile_idx,
                log_sum!(
                    dp_matrix.get_match(target_idx + 1, profile_idx + 1)
                        + profile.transition_score(Profile::DELETE_TO_MATCH_IDX, profile_idx)
                        + profile.match_score(current_residue, profile_idx + 1),
                    dp_matrix.get_delete(target_idx, profile_idx + 1)
                        + profile.transition_score(Profile::DELETE_TO_DELETE_IDX, profile_idx),
                    dp_matrix.get_special(target_idx, Profile::SPECIAL_E_IDX) + end_score
                ),
            );
        }
    }

    let first_target_character = target.digital_bytes[row_bounds.target_start] as usize;

    let profile_start_in_first_row = row_bounds.left_row_bounds[row_bounds.target_start];
    let profile_end_in_first_row = row_bounds.right_row_bounds[row_bounds.target_start];

    dp_matrix.set_special(
        row_bounds.target_start - 1,
        Profile::SPECIAL_B_IDX,
        dp_matrix.get_match(row_bounds.target_start, profile_start_in_first_row)
            + profile.transition_score(Profile::BEGIN_TO_MATCH_IDX, 0)
            + profile.match_score(first_target_character, 1),
    );

    for profile_idx in (profile_start_in_first_row + 1)..=profile_end_in_first_row {
        dp_matrix.set_special(
            row_bounds.target_start - 1,
            Profile::SPECIAL_B_IDX,
            log_sum!(
                dp_matrix.get_special(row_bounds.target_start - 1, Profile::SPECIAL_B_IDX),
                dp_matrix.get_match(row_bounds.target_start, profile_idx)
                    + profile.transition_score(Profile::BEGIN_TO_MATCH_IDX, profile_idx - 1)
                    + profile.match_score(first_target_character, profile_idx)
            ),
        );
    }

    dp_matrix.set_special(
        row_bounds.target_start - 1,
        Profile::SPECIAL_J_IDX,
        -f32::INFINITY,
    );
    dp_matrix.set_special(
        row_bounds.target_start - 1,
        Profile::SPECIAL_C_IDX,
        -f32::INFINITY,
    );
    dp_matrix.set_special(
        row_bounds.target_start - 1,
        Profile::SPECIAL_E_IDX,
        -f32::INFINITY,
    );
    dp_matrix.set_special(
        row_bounds.target_start - 1,
        Profile::SPECIAL_N_IDX,
        log_sum!(
            dp_matrix.get_special(row_bounds.target_start, Profile::SPECIAL_N_IDX)
                + profile
                    .special_transition_score(Profile::SPECIAL_N_IDX, Profile::SPECIAL_LOOP_IDX),
            dp_matrix.get_special(row_bounds.target_start - 1, Profile::SPECIAL_B_IDX)
                + profile
                    .special_transition_score(Profile::SPECIAL_N_IDX, Profile::SPECIAL_MOVE_IDX)
        ),
    );
    for profile_idx in (profile_start_in_first_row..=profile_end_in_first_row).rev() {
        dp_matrix.set_match(row_bounds.target_start - 1, profile_idx, -f32::INFINITY);
        dp_matrix.set_insert(row_bounds.target_start - 1, profile_idx, -f32::INFINITY);
        dp_matrix.set_delete(row_bounds.target_start - 1, profile_idx, -f32::INFINITY);
    }
}
