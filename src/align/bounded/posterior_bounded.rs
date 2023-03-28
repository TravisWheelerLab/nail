use crate::align::bounded::structs::RowBoundParams;
use crate::structs::dp_matrix::DpMatrix;
use crate::structs::Profile;
use crate::timing::time;

#[funci::timed(timer = time)]
pub fn posterior_bounded(
    profile: &Profile,
    // forward_matrix: &DpMatrix3D,
    // backward_matrix: &DpMatrix3D,
    // posterior_matrix: &mut DpMatrix3D,
    forward_matrix: &impl DpMatrix,
    backward_matrix: &impl DpMatrix,
    posterior_matrix: &mut impl DpMatrix,
    params: &RowBoundParams,
) {
    // let target_length = forward_matrix.target_length;

    let overall_score: f32 = forward_matrix.get_special(params.target_end, Profile::SPECIAL_C_IDX)
        + profile.special_transition_score(Profile::SPECIAL_C_IDX, Profile::SPECIAL_MOVE_IDX);

    let mut denominator: f32;

    posterior_matrix.set_special(params.target_start - 1, Profile::SPECIAL_E_IDX, 0.0);
    posterior_matrix.set_special(params.target_start - 1, Profile::SPECIAL_N_IDX, 0.0);
    posterior_matrix.set_special(params.target_start - 1, Profile::SPECIAL_J_IDX, 0.0);
    posterior_matrix.set_special(params.target_start - 1, Profile::SPECIAL_B_IDX, 0.0);
    posterior_matrix.set_special(params.target_start - 1, Profile::SPECIAL_C_IDX, 0.0);

    let profile_start_in_first_row = params.left_row_bounds[params.target_start];
    let profile_end_in_first_row = params.right_row_bounds[params.target_start];

    // for profile_idx in 0..=profile.length {
    for profile_idx in (profile_start_in_first_row - 1)..=profile_end_in_first_row {
        posterior_matrix.set_match(params.target_start - 1, profile_idx, 0.0);
        posterior_matrix.set_insert(params.target_start - 1, profile_idx, 0.0);
        posterior_matrix.set_delete(params.target_start - 1, profile_idx, 0.0);
    }

    // for target_idx in 1..=target_length {
    for target_idx in params.target_start..=params.target_end {
        denominator = 0.0;
        posterior_matrix.set_match(target_idx, params.target_start - 1, 0.0);
        posterior_matrix.set_insert(target_idx, params.target_start - 1, 0.0);
        posterior_matrix.set_delete(target_idx, params.target_start - 1, 0.0);

        let profile_start_in_current_row = params.left_row_bounds[target_idx];
        let profile_end_in_current_row = params.right_row_bounds[target_idx];

        // for profile_idx in 1..profile.length {
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

        // for profile_idx in 1..profile.length {
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
