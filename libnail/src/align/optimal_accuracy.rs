use crate::align::structs::{DpMatrix, RowBounds};
use crate::max_f32;
use crate::structs::Profile;

pub fn optimal_accuracy(
    profile: &Profile,
    posterior_matrix: &impl DpMatrix,
    optimal_matrix: &mut impl DpMatrix,
    row_bounds: &RowBounds,
) {
    let end_score: f32 = 1.0;

    // initialization of the zero row
    optimal_matrix.set_special(row_bounds.target_start - 1, Profile::SPECIAL_N_IDX, 0.0);
    optimal_matrix.set_special(row_bounds.target_start - 1, Profile::SPECIAL_B_IDX, 0.0);
    optimal_matrix.set_special(
        row_bounds.target_start - 1,
        Profile::SPECIAL_E_IDX,
        -f32::INFINITY,
    );
    optimal_matrix.set_special(
        row_bounds.target_start - 1,
        Profile::SPECIAL_C_IDX,
        -f32::INFINITY,
    );
    optimal_matrix.set_special(
        row_bounds.target_start - 1,
        Profile::SPECIAL_J_IDX,
        -f32::INFINITY,
    );

    let profile_start_in_first_row = row_bounds.left_row_bounds[row_bounds.target_start];
    let profile_end_in_first_row = row_bounds.right_row_bounds[row_bounds.target_start];

    // for profile_idx in 0..=profile.length {
    for profile_idx in (profile_start_in_first_row - 1)..=profile_end_in_first_row {
        optimal_matrix.set_match(row_bounds.target_start - 1, profile_idx, -f32::INFINITY);
        optimal_matrix.set_insert(row_bounds.target_start - 1, profile_idx, -f32::INFINITY);
        optimal_matrix.set_delete(row_bounds.target_start - 1, profile_idx, -f32::INFINITY);
    }

    // for i in 1..=posterior_matrix.target_length {
    for target_idx in row_bounds.target_start..=row_bounds.target_end {
        let profile_start_in_current_row = row_bounds.left_row_bounds[target_idx];
        let profile_end_in_current_row = row_bounds.right_row_bounds[target_idx];

        optimal_matrix.set_match(target_idx, profile_start_in_current_row - 1, -f32::INFINITY);
        optimal_matrix.set_insert(target_idx, profile_start_in_current_row - 1, -f32::INFINITY);
        optimal_matrix.set_delete(target_idx, profile_start_in_current_row - 1, -f32::INFINITY);
        optimal_matrix.set_special(target_idx, Profile::SPECIAL_E_IDX, -f32::INFINITY);

        // for profile_idx in 1..profile.length {
        for profile_idx in profile_start_in_current_row..profile_end_in_current_row {
            optimal_matrix.set_match(
                target_idx,
                profile_idx,
                max_f32!(
                    profile.transition_score_delta(Profile::M_M_IDX, profile_idx - 1)
                        * (optimal_matrix.get_match(target_idx - 1, profile_idx - 1)
                            + posterior_matrix.get_match(target_idx, profile_idx)),
                    profile.transition_score_delta(Profile::I_M_IDX, profile_idx - 1)
                        * (optimal_matrix.get_insert(target_idx - 1, profile_idx - 1)
                            + posterior_matrix.get_match(target_idx, profile_idx)),
                    profile.transition_score_delta(Profile::D_M_IDX, profile_idx - 1)
                        * (optimal_matrix.get_delete(target_idx - 1, profile_idx - 1)
                            + posterior_matrix.get_match(target_idx, profile_idx)),
                    profile.transition_score_delta(Profile::B_M_IDX, profile_idx - 1)
                        * (optimal_matrix.get_special(target_idx - 1, Profile::SPECIAL_B_IDX)
                            + posterior_matrix.get_match(target_idx, profile_idx))
                ),
            );

            optimal_matrix.set_special(
                target_idx,
                Profile::SPECIAL_E_IDX,
                max_f32!(
                    optimal_matrix.get_special(target_idx, Profile::SPECIAL_E_IDX),
                    optimal_matrix.get_match(target_idx, profile_idx) * end_score
                ),
            );

            optimal_matrix.set_insert(
                target_idx,
                profile_idx,
                max_f32!(
                    profile.transition_score_delta(Profile::M_I_IDX, profile_idx)
                        * (optimal_matrix.get_match(target_idx - 1, profile_idx)
                            + posterior_matrix.get_insert(target_idx, profile_idx)),
                    profile.transition_score_delta(Profile::I_I_IDX, profile_idx)
                        * (optimal_matrix.get_insert(target_idx - 1, profile_idx)
                            + posterior_matrix.get_insert(target_idx, profile_idx))
                ),
            );

            optimal_matrix.set_delete(
                target_idx,
                profile_idx,
                max_f32!(
                    profile.transition_score_delta(Profile::M_D_IDX, profile_idx - 1)
                        * optimal_matrix.get_match(target_idx, profile_idx - 1),
                    profile.transition_score_delta(Profile::D_D_IDX, profile_idx - 1)
                        * optimal_matrix.get_delete(target_idx, profile_idx - 1)
                ),
            );
        }

        optimal_matrix.set_match(
            target_idx,
            profile_end_in_current_row,
            max_f32!(
                profile.transition_score_delta(Profile::M_M_IDX, profile_end_in_current_row - 1)
                    * (optimal_matrix.get_match(target_idx - 1, profile_end_in_current_row - 1)
                        + posterior_matrix.get_match(target_idx, profile_end_in_current_row)),
                profile.transition_score_delta(Profile::I_M_IDX, profile_end_in_current_row - 1)
                    * (optimal_matrix.get_insert(target_idx - 1, profile_end_in_current_row - 1)
                        + posterior_matrix.get_match(target_idx, profile_end_in_current_row)),
                profile.transition_score_delta(Profile::D_M_IDX, profile_end_in_current_row - 1)
                    * (optimal_matrix.get_delete(target_idx - 1, profile_end_in_current_row - 1)
                        + posterior_matrix.get_match(target_idx, profile_end_in_current_row)),
                profile.transition_score_delta(Profile::B_M_IDX, profile_end_in_current_row - 1)
                    * (optimal_matrix.get_special(target_idx - 1, Profile::SPECIAL_B_IDX)
                        + posterior_matrix.get_match(target_idx, profile_end_in_current_row))
            ),
        );

        optimal_matrix.set_delete(
            target_idx,
            profile_end_in_current_row,
            max_f32!(
                profile.transition_score_delta(Profile::M_D_IDX, profile_end_in_current_row - 1)
                    * optimal_matrix.get_match(target_idx, profile_end_in_current_row - 1),
                profile.transition_score_delta(Profile::D_D_IDX, profile_end_in_current_row - 1)
                    * optimal_matrix.get_delete(target_idx, profile_end_in_current_row - 1)
            ),
        );

        // a comment from hmmer:
        //   now the special states; it's important that E is already done, and B is done after N,J
        optimal_matrix.set_special(
            target_idx,
            Profile::SPECIAL_E_IDX,
            max_f32!(
                optimal_matrix.get_special(target_idx, Profile::SPECIAL_E_IDX),
                optimal_matrix.get_match(target_idx, profile_end_in_current_row),
                optimal_matrix.get_delete(target_idx, profile_end_in_current_row)
            ),
        );

        optimal_matrix.set_special(
            target_idx,
            Profile::SPECIAL_J_IDX,
            max_f32!(
                profile.special_transition_score_delta(
                    Profile::SPECIAL_J_IDX,
                    Profile::SPECIAL_LOOP_IDX
                ) * (optimal_matrix.get_special(target_idx - 1, Profile::SPECIAL_J_IDX)
                    + posterior_matrix.get_special(target_idx, Profile::SPECIAL_J_IDX)),
                profile.special_transition_score_delta(
                    Profile::SPECIAL_E_IDX,
                    Profile::SPECIAL_LOOP_IDX
                ) * optimal_matrix.get_special(target_idx, Profile::SPECIAL_E_IDX)
            ),
        );

        optimal_matrix.set_special(
            target_idx,
            Profile::SPECIAL_C_IDX,
            max_f32!(
                profile.special_transition_score_delta(
                    Profile::SPECIAL_C_IDX,
                    Profile::SPECIAL_LOOP_IDX
                ) * (optimal_matrix.get_special(target_idx - 1, Profile::SPECIAL_C_IDX)
                    + posterior_matrix.get_special(target_idx, Profile::SPECIAL_C_IDX)),
                profile.special_transition_score_delta(
                    Profile::SPECIAL_E_IDX,
                    Profile::SPECIAL_MOVE_IDX
                ) * optimal_matrix.get_special(target_idx, Profile::SPECIAL_E_IDX)
            ),
        );

        optimal_matrix.set_special(
            target_idx,
            Profile::SPECIAL_N_IDX,
            profile
                .special_transition_score_delta(Profile::SPECIAL_N_IDX, Profile::SPECIAL_LOOP_IDX)
                * (optimal_matrix.get_special(target_idx - 1, Profile::SPECIAL_N_IDX)
                    + posterior_matrix.get_special(target_idx, Profile::SPECIAL_N_IDX)),
        );

        optimal_matrix.set_special(
            target_idx,
            Profile::SPECIAL_B_IDX,
            max_f32!(
                profile.special_transition_score_delta(
                    Profile::SPECIAL_N_IDX,
                    Profile::SPECIAL_MOVE_IDX
                ) * optimal_matrix.get_special(target_idx, Profile::SPECIAL_N_IDX),
                profile.special_transition_score_delta(
                    Profile::SPECIAL_J_IDX,
                    Profile::SPECIAL_MOVE_IDX
                ) * optimal_matrix.get_special(target_idx, Profile::SPECIAL_J_IDX)
            ),
        );
    }
}
