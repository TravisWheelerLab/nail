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
    optimal_matrix.set_special(0, Profile::SPECIAL_N, 0.0);
    optimal_matrix.set_special(0, Profile::SPECIAL_B, 0.0);
    optimal_matrix.set_special(0, Profile::SPECIAL_E, -f32::INFINITY);
    optimal_matrix.set_special(0, Profile::SPECIAL_C, -f32::INFINITY);
    optimal_matrix.set_special(0, Profile::SPECIAL_J, -f32::INFINITY);

    for profile_idx in 0..=profile.length {
        optimal_matrix.set_match(0, profile_idx, -f32::INFINITY);
        optimal_matrix.set_insert(0, profile_idx, -f32::INFINITY);
        optimal_matrix.set_delete(0, profile_idx, -f32::INFINITY);
    }

    for i in 1..=posterior_matrix.target_length() {
        optimal_matrix.set_match(i, 0, -f32::INFINITY);
        optimal_matrix.set_insert(i, 0, -f32::INFINITY);
        optimal_matrix.set_delete(i, 0, -f32::INFINITY);
        optimal_matrix.set_special(i, Profile::SPECIAL_E, -f32::INFINITY);

        for k in 1..profile.length {
            optimal_matrix.set_match(
                i,
                k,
                f32::max(
                    f32::max(
                        profile.transition_score_delta(Profile::PROFILE_MATCH_TO_MATCH, k - 1)
                            * (optimal_matrix.get_match(i - 1, k - 1)
                                + posterior_matrix.get_match(i, k)),
                        profile.transition_score_delta(Profile::PROFILE_INSERT_TO_MATCH, k - 1)
                            * (optimal_matrix.get_insert(i - 1, k - 1)
                                + posterior_matrix.get_match(i, k)),
                    ),
                    f32::max(
                        profile.transition_score_delta(Profile::PROFILE_DELETE_TO_MATCH, k - 1)
                            * (optimal_matrix.get_delete(i - 1, k - 1)
                                + posterior_matrix.get_match(i, k)),
                        profile.transition_score_delta(Profile::PROFILE_BEGIN_TO_MATCH, k - 1)
                            * (optimal_matrix.get_special(i - 1, Profile::SPECIAL_B)
                                + posterior_matrix.get_match(i, k)),
                    ),
                ),
            );

            // E computed before delete here as an alleged optimization
            // because we probably had the match_score(i, k) on a register
            // TODO: see if this actually makes any measurable difference
            optimal_matrix.set_special(
                i,
                Profile::SPECIAL_E,
                f32::max(
                    optimal_matrix.get_special(i, Profile::SPECIAL_E),
                    optimal_matrix.get_match(i, k) * end_score,
                ),
            );

            optimal_matrix.set_insert(
                i,
                k,
                f32::max(
                    profile.transition_score_delta(Profile::PROFILE_MATCH_TO_INSERT, k)
                        * (optimal_matrix.get_match(i - 1, k) + posterior_matrix.get_insert(i, k)),
                    profile.transition_score_delta(Profile::PROFILE_INSERT_TO_INSERT, k)
                        * (optimal_matrix.get_insert(i - 1, k) + posterior_matrix.get_insert(i, k)),
                ),
            );

            optimal_matrix.set_delete(
                i,
                k,
                f32::max(
                    profile.transition_score_delta(Profile::PROFILE_MATCH_TO_DELETE, k - 1)
                        * optimal_matrix.get_match(i, k - 1),
                    profile.transition_score_delta(Profile::PROFILE_DELETE_TO_DELETE, k - 1)
                        * optimal_matrix.get_delete(i, k - 1),
                ),
            );
        }

        optimal_matrix.set_match(
            i,
            profile.length,
            f32::max(
                f32::max(
                    profile.transition_score_delta(
                        Profile::PROFILE_MATCH_TO_MATCH,
                        profile.length - 1,
                    ) * (optimal_matrix.get_match(i - 1, profile.length - 1)
                        + posterior_matrix.get_match(i, profile.length)),
                    profile.transition_score_delta(
                        Profile::PROFILE_INSERT_TO_MATCH,
                        profile.length - 1,
                    ) * (optimal_matrix.get_insert(i - 1, profile.length - 1)
                        + posterior_matrix.get_match(i, profile.length)),
                ),
                f32::max(
                    profile.transition_score_delta(
                        Profile::PROFILE_DELETE_TO_MATCH,
                        profile.length - 1,
                    ) * (optimal_matrix.get_delete(i - 1, profile.length - 1)
                        + posterior_matrix.get_match(i, profile.length)),
                    profile.transition_score_delta(
                        Profile::PROFILE_BEGIN_TO_MATCH,
                        profile.length - 1,
                    ) * (optimal_matrix.get_special(i - 1, Profile::SPECIAL_B)
                        + posterior_matrix.get_match(i, profile.length)),
                ),
            ),
        );

        optimal_matrix.set_delete(
            i,
            profile.length,
            f32::max(
                profile
                    .transition_score_delta(Profile::PROFILE_MATCH_TO_DELETE, profile.length - 1)
                    * optimal_matrix.get_match(i, profile.length - 1),
                profile
                    .transition_score_delta(Profile::PROFILE_DELETE_TO_DELETE, profile.length - 1)
                    * optimal_matrix.get_delete(i, profile.length - 1),
            ),
        );

        // a comment from hmmer:
        //   now the special states; it's important that E is already done, and B is done after N,J
        optimal_matrix.set_special(
            i,
            Profile::SPECIAL_E,
            f32::max(
                optimal_matrix.get_special(i, Profile::SPECIAL_E),
                f32::max(
                    optimal_matrix.get_match(i, profile.length),
                    optimal_matrix.get_delete(i, profile.length),
                ),
            ),
        );

        optimal_matrix.set_special(
            i,
            Profile::SPECIAL_J,
            f32::max(
                profile.special_transition_score_delta(Profile::SPECIAL_J, Profile::SPECIAL_LOOP)
                    * (optimal_matrix.get_special(i - 1, Profile::SPECIAL_J)
                        + posterior_matrix.get_special(i, Profile::SPECIAL_J)),
                profile.special_transition_score_delta(Profile::SPECIAL_E, Profile::SPECIAL_LOOP)
                    * optimal_matrix.get_special(i, Profile::SPECIAL_E),
            ),
        );

        optimal_matrix.set_special(
            i,
            Profile::SPECIAL_C,
            f32::max(
                profile.special_transition_score_delta(Profile::SPECIAL_C, Profile::SPECIAL_LOOP)
                    * (optimal_matrix.get_special(i - 1, Profile::SPECIAL_C)
                        + posterior_matrix.get_special(i, Profile::SPECIAL_C)),
                profile.special_transition_score_delta(Profile::SPECIAL_E, Profile::SPECIAL_MOVE)
                    * optimal_matrix.get_special(i, Profile::SPECIAL_E),
            ),
        );

        optimal_matrix.set_special(
            i,
            Profile::SPECIAL_N,
            profile.special_transition_score_delta(Profile::SPECIAL_N, Profile::SPECIAL_LOOP)
                * (optimal_matrix.get_special(i - 1, Profile::SPECIAL_N)
                    + posterior_matrix.get_special(i, Profile::SPECIAL_N)),
        );

        optimal_matrix.set_special(
            i,
            Profile::SPECIAL_B,
            f32::max(
                profile.special_transition_score_delta(Profile::SPECIAL_N, Profile::SPECIAL_MOVE)
                    * optimal_matrix.get_special(i, Profile::SPECIAL_N),
                profile.special_transition_score_delta(Profile::SPECIAL_J, Profile::SPECIAL_MOVE)
                    * optimal_matrix.get_special(i, Profile::SPECIAL_J),
            ),
        );
    }
}
