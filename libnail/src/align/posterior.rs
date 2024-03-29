use crate::align::structs::{DpMatrix, RowBounds};
use crate::structs::Profile;

pub fn posterior(
    profile: &Profile,
    forward_matrix: &impl DpMatrix,
    backward_matrix: &impl DpMatrix,
    posterior_matrix: &mut impl DpMatrix,
    row_bounds: &RowBounds,
) {
    let overall_score: f32 = forward_matrix
        .get_special(row_bounds.target_end, Profile::SPECIAL_C_IDX)
        + profile.special_transition_score(Profile::SPECIAL_C_IDX, Profile::SPECIAL_MOVE_IDX);

    let mut denominator: f32;

    posterior_matrix.set_special(row_bounds.target_start - 1, Profile::SPECIAL_E_IDX, 0.0);
    posterior_matrix.set_special(row_bounds.target_start - 1, Profile::SPECIAL_N_IDX, 0.0);
    posterior_matrix.set_special(row_bounds.target_start - 1, Profile::SPECIAL_J_IDX, 0.0);
    posterior_matrix.set_special(row_bounds.target_start - 1, Profile::SPECIAL_B_IDX, 0.0);
    posterior_matrix.set_special(row_bounds.target_start - 1, Profile::SPECIAL_C_IDX, 0.0);

    let profile_start_in_first_row = row_bounds.left_row_bounds[row_bounds.target_start];
    let profile_end_in_first_row = row_bounds.right_row_bounds[row_bounds.target_start];

    for profile_idx in (profile_start_in_first_row - 1)..=profile_end_in_first_row {
        posterior_matrix.set_match(row_bounds.target_start - 1, profile_idx, 0.0);
        posterior_matrix.set_insert(row_bounds.target_start - 1, profile_idx, 0.0);
        posterior_matrix.set_delete(row_bounds.target_start - 1, profile_idx, 0.0);
    }

    for target_idx in row_bounds.target_start..=row_bounds.target_end {
        denominator = 0.0;

        let profile_start_in_current_row = row_bounds.left_row_bounds[target_idx];
        let profile_end_in_current_row = row_bounds.right_row_bounds[target_idx];

        posterior_matrix.set_match(target_idx, profile_start_in_current_row - 1, 0.0);
        posterior_matrix.set_insert(target_idx, profile_start_in_current_row - 1, 0.0);
        posterior_matrix.set_delete(target_idx, profile_start_in_current_row - 1, 0.0);

        for profile_idx in profile_start_in_current_row..profile_end_in_current_row {
            posterior_matrix.set_match(
                target_idx,
                profile_idx,
                (forward_matrix.get_match(target_idx, profile_idx)
                    + backward_matrix.get_match(target_idx, profile_idx)
                    - overall_score)
                    .exp(),
            );

            denominator += posterior_matrix.get_match(target_idx, profile_idx);

            posterior_matrix.set_insert(
                target_idx,
                profile_idx,
                (forward_matrix.get_insert(target_idx, profile_idx)
                    + backward_matrix.get_insert(target_idx, profile_idx)
                    - overall_score)
                    .exp(),
            );

            denominator += posterior_matrix.get_insert(target_idx, profile_idx);

            posterior_matrix.set_delete(target_idx, profile_idx, 0.0);
        }

        posterior_matrix.set_match(
            target_idx,
            profile_end_in_current_row,
            (forward_matrix.get_match(target_idx, profile_end_in_current_row)
                + backward_matrix.get_match(target_idx, profile_end_in_current_row)
                - overall_score)
                .exp(),
        );

        denominator += posterior_matrix.get_match(target_idx, profile_end_in_current_row);
        posterior_matrix.set_insert(target_idx, profile_end_in_current_row, 0.0);
        posterior_matrix.set_delete(target_idx, profile_end_in_current_row, 0.0);

        posterior_matrix.set_special(target_idx, Profile::SPECIAL_E_IDX, 0.0);

        posterior_matrix.set_special(
            target_idx,
            Profile::SPECIAL_N_IDX,
            (forward_matrix.get_special(target_idx - 1, Profile::SPECIAL_N_IDX)
                + backward_matrix.get_special(target_idx, Profile::SPECIAL_N_IDX)
                + profile
                    .special_transition_score(Profile::SPECIAL_N_IDX, Profile::SPECIAL_LOOP_IDX)
                - overall_score)
                .exp(),
        );

        posterior_matrix.set_special(
            target_idx,
            Profile::SPECIAL_J_IDX,
            (forward_matrix.get_special(target_idx - 1, Profile::SPECIAL_J_IDX)
                + backward_matrix.get_special(target_idx, Profile::SPECIAL_J_IDX)
                + profile
                    .special_transition_score(Profile::SPECIAL_J_IDX, Profile::SPECIAL_LOOP_IDX)
                - overall_score)
                .exp(),
        );

        posterior_matrix.set_special(target_idx, Profile::SPECIAL_B_IDX, 0.0);

        posterior_matrix.set_special(
            target_idx,
            Profile::SPECIAL_C_IDX,
            (forward_matrix.get_special(target_idx - 1, Profile::SPECIAL_C_IDX)
                + backward_matrix.get_special(target_idx, Profile::SPECIAL_C_IDX)
                + profile
                    .special_transition_score(Profile::SPECIAL_C_IDX, Profile::SPECIAL_LOOP_IDX)
                - overall_score)
                .exp(),
        );

        denominator += posterior_matrix.get_special(target_idx, Profile::SPECIAL_N_IDX);
        denominator += posterior_matrix.get_special(target_idx, Profile::SPECIAL_J_IDX);
        denominator += posterior_matrix.get_special(target_idx, Profile::SPECIAL_C_IDX);

        denominator = 1.0 / denominator;

        for profile_idx in profile_start_in_current_row..profile_end_in_current_row {
            posterior_matrix.set_match(
                target_idx,
                profile_idx,
                posterior_matrix.get_match(target_idx, profile_idx) * denominator,
            );

            posterior_matrix.set_insert(
                target_idx,
                profile_idx,
                posterior_matrix.get_insert(target_idx, profile_idx) * denominator,
            );
        }

        posterior_matrix.set_match(
            target_idx,
            profile_end_in_current_row,
            posterior_matrix.get_match(target_idx, profile_end_in_current_row) * denominator,
        );

        posterior_matrix.set_special(
            target_idx,
            Profile::SPECIAL_N_IDX,
            posterior_matrix.get_special(target_idx, Profile::SPECIAL_N_IDX) * denominator,
        );

        posterior_matrix.set_special(
            target_idx,
            Profile::SPECIAL_J_IDX,
            posterior_matrix.get_special(target_idx, Profile::SPECIAL_J_IDX) * denominator,
        );

        posterior_matrix.set_special(
            target_idx,
            Profile::SPECIAL_C_IDX,
            posterior_matrix.get_special(target_idx, Profile::SPECIAL_C_IDX) * denominator,
        );
    }
}
