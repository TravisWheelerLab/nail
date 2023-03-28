use crate::structs::dp_matrix::DpMatrix;
use crate::structs::Profile;
use crate::timing::time;

#[funci::timed(timer = time)]
pub fn optimal_accuracy(
    profile: &Profile,
    posterior_matrix: &impl DpMatrix,
    optimal_matrix: &mut impl DpMatrix,
) {
    let end_score: f32 = 1.0;

    // initialization of the zero row
    optimal_matrix.set_special(0, Profile::SPECIAL_N_IDX, 0.0);
    optimal_matrix.set_special(0, Profile::SPECIAL_B_IDX, 0.0);
    optimal_matrix.set_special(0, Profile::SPECIAL_E_IDX, -f32::INFINITY);
    optimal_matrix.set_special(0, Profile::SPECIAL_C_IDX, -f32::INFINITY);
    optimal_matrix.set_special(0, Profile::SPECIAL_J_IDX, -f32::INFINITY);

    for profile_idx in 0..=profile.length {
        optimal_matrix.set_match(0, profile_idx, -f32::INFINITY);
        optimal_matrix.set_insert(0, profile_idx, -f32::INFINITY);
        optimal_matrix.set_delete(0, profile_idx, -f32::INFINITY);
    }

    for i in 1..=posterior_matrix.target_length() {
        optimal_matrix.set_match(i, 0, -f32::INFINITY);
        optimal_matrix.set_insert(i, 0, -f32::INFINITY);
        optimal_matrix.set_delete(i, 0, -f32::INFINITY);
        optimal_matrix.set_special(i, Profile::SPECIAL_E_IDX, -f32::INFINITY);

        for k in 1..profile.length {
            optimal_matrix.set_match(
                i,
                k,
                f32::max(
                    f32::max(
                        profile.transition_score_delta(Profile::MATCH_TO_MATCH_IDX, k - 1)
                            * (optimal_matrix.get_match(i - 1, k - 1)
                                + posterior_matrix.get_match(i, k)),
                        profile.transition_score_delta(Profile::INSERT_TO_MATCH_IDX, k - 1)
                            * (optimal_matrix.get_insert(i - 1, k - 1)
                                + posterior_matrix.get_match(i, k)),
                    ),
                    f32::max(
                        profile.transition_score_delta(Profile::DELETE_TO_MATCH_IDX, k - 1)
                            * (optimal_matrix.get_delete(i - 1, k - 1)
                                + posterior_matrix.get_match(i, k)),
                        profile.transition_score_delta(Profile::BEGIN_TO_MATCH_IDX, k - 1)
                            * (optimal_matrix.get_special(i - 1, Profile::SPECIAL_B_IDX)
                                + posterior_matrix.get_match(i, k)),
                    ),
                ),
            );

            // E computed before delete here as an alleged optimization
            // because we probably had the match_score(i, k) on a register
            // TODO: see if this actually makes any measurable difference
            optimal_matrix.set_special(
                i,
                Profile::SPECIAL_E_IDX,
                f32::max(
                    optimal_matrix.get_special(i, Profile::SPECIAL_E_IDX),
                    optimal_matrix.get_match(i, k) * end_score,
                ),
            );

            optimal_matrix.set_insert(
                i,
                k,
                f32::max(
                    profile.transition_score_delta(Profile::MATCH_TO_INSERT_IDX, k)
                        * (optimal_matrix.get_match(i - 1, k) + posterior_matrix.get_insert(i, k)),
                    profile.transition_score_delta(Profile::INSERT_TO_INSERT_IDX, k)
                        * (optimal_matrix.get_insert(i - 1, k) + posterior_matrix.get_insert(i, k)),
                ),
            );

            optimal_matrix.set_delete(
                i,
                k,
                f32::max(
                    profile.transition_score_delta(Profile::MATCH_TO_DELETE_IDX, k - 1)
                        * optimal_matrix.get_match(i, k - 1),
                    profile.transition_score_delta(Profile::DELETE_TO_DELETE_IDX, k - 1)
                        * optimal_matrix.get_delete(i, k - 1),
                ),
            );
        }

        optimal_matrix.set_match(
            i,
            profile.length,
            f32::max(
                f32::max(
                    profile.transition_score_delta(Profile::MATCH_TO_MATCH_IDX, profile.length - 1)
                        * (optimal_matrix.get_match(i - 1, profile.length - 1)
                            + posterior_matrix.get_match(i, profile.length)),
                    profile
                        .transition_score_delta(Profile::INSERT_TO_MATCH_IDX, profile.length - 1)
                        * (optimal_matrix.get_insert(i - 1, profile.length - 1)
                            + posterior_matrix.get_match(i, profile.length)),
                ),
                f32::max(
                    profile
                        .transition_score_delta(Profile::DELETE_TO_MATCH_IDX, profile.length - 1)
                        * (optimal_matrix.get_delete(i - 1, profile.length - 1)
                            + posterior_matrix.get_match(i, profile.length)),
                    profile.transition_score_delta(Profile::BEGIN_TO_MATCH_IDX, profile.length - 1)
                        * (optimal_matrix.get_special(i - 1, Profile::SPECIAL_B_IDX)
                            + posterior_matrix.get_match(i, profile.length)),
                ),
            ),
        );

        optimal_matrix.set_delete(
            i,
            profile.length,
            f32::max(
                profile.transition_score_delta(Profile::MATCH_TO_DELETE_IDX, profile.length - 1)
                    * optimal_matrix.get_match(i, profile.length - 1),
                profile.transition_score_delta(Profile::DELETE_TO_DELETE_IDX, profile.length - 1)
                    * optimal_matrix.get_delete(i, profile.length - 1),
            ),
        );

        // a comment from hmmer:
        //   now the special states; it's important that E is already done, and B is done after N,J
        optimal_matrix.set_special(
            i,
            Profile::SPECIAL_E_IDX,
            f32::max(
                optimal_matrix.get_special(i, Profile::SPECIAL_E_IDX),
                f32::max(
                    optimal_matrix.get_match(i, profile.length),
                    optimal_matrix.get_delete(i, profile.length),
                ),
            ),
        );

        optimal_matrix.set_special(
            i,
            Profile::SPECIAL_J_IDX,
            f32::max(
                profile.special_transition_score_delta(
                    Profile::SPECIAL_J_IDX,
                    Profile::SPECIAL_LOOP_IDX,
                ) * (optimal_matrix.get_special(i - 1, Profile::SPECIAL_J_IDX)
                    + posterior_matrix.get_special(i, Profile::SPECIAL_J_IDX)),
                profile.special_transition_score_delta(
                    Profile::SPECIAL_E_IDX,
                    Profile::SPECIAL_LOOP_IDX,
                ) * optimal_matrix.get_special(i, Profile::SPECIAL_E_IDX),
            ),
        );

        optimal_matrix.set_special(
            i,
            Profile::SPECIAL_C_IDX,
            f32::max(
                profile.special_transition_score_delta(
                    Profile::SPECIAL_C_IDX,
                    Profile::SPECIAL_LOOP_IDX,
                ) * (optimal_matrix.get_special(i - 1, Profile::SPECIAL_C_IDX)
                    + posterior_matrix.get_special(i, Profile::SPECIAL_C_IDX)),
                profile.special_transition_score_delta(
                    Profile::SPECIAL_E_IDX,
                    Profile::SPECIAL_MOVE_IDX,
                ) * optimal_matrix.get_special(i, Profile::SPECIAL_E_IDX),
            ),
        );

        optimal_matrix.set_special(
            i,
            Profile::SPECIAL_N_IDX,
            profile
                .special_transition_score_delta(Profile::SPECIAL_N_IDX, Profile::SPECIAL_LOOP_IDX)
                * (optimal_matrix.get_special(i - 1, Profile::SPECIAL_N_IDX)
                    + posterior_matrix.get_special(i, Profile::SPECIAL_N_IDX)),
        );

        optimal_matrix.set_special(
            i,
            Profile::SPECIAL_B_IDX,
            f32::max(
                profile.special_transition_score_delta(
                    Profile::SPECIAL_N_IDX,
                    Profile::SPECIAL_MOVE_IDX,
                ) * optimal_matrix.get_special(i, Profile::SPECIAL_N_IDX),
                profile.special_transition_score_delta(
                    Profile::SPECIAL_J_IDX,
                    Profile::SPECIAL_MOVE_IDX,
                ) * optimal_matrix.get_special(i, Profile::SPECIAL_J_IDX),
            ),
        );
    }
}
