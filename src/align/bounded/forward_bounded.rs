use crate::align::bounded::structs::RowBoundParams;
use crate::log_sum;
use crate::structs::dp_matrix::DpMatrix;
use crate::structs::{Profile, Sequence};
use crate::timing::time;
use crate::util::log_add;

#[funci::timed(timer = time)]
pub fn forward_bounded(
    profile: &Profile,
    target: &Sequence,
    dp_matrix: &mut impl DpMatrix,
    params: &RowBoundParams,
) {
    let end_score: f32 = 0.0;

    dp_matrix.set_special(params.target_start - 1, Profile::SPECIAL_N_IDX, 0.0);
    dp_matrix.set_special(
        params.target_start - 1,
        Profile::SPECIAL_B_IDX,
        profile.special_transition_score(Profile::SPECIAL_N_IDX, Profile::SPECIAL_MOVE_IDX),
    );
    dp_matrix.set_special(
        params.target_start - 1,
        Profile::SPECIAL_E_IDX,
        -f32::INFINITY,
    );
    dp_matrix.set_special(
        params.target_start - 1,
        Profile::SPECIAL_C_IDX,
        -f32::INFINITY,
    );
    dp_matrix.set_special(
        params.target_start - 1,
        Profile::SPECIAL_J_IDX,
        -f32::INFINITY,
    );

    // TODO: probably need B & N state computation from 0..target_start
    //       also maybe remove J state?

    for target_idx in params.target_start..=params.target_end {
        let current_target_character = target.digital_bytes[target_idx];

        for profile_idx in params.left_row_bounds[target_idx]..params.right_row_bounds[target_idx] {
            // match state
            dp_matrix.set_match(
                target_idx,
                profile_idx,
                log_sum!(
                    dp_matrix.get_match(target_idx - 1, profile_idx - 1)
                        + profile.transition_score(Profile::MATCH_TO_MATCH_IDX, profile_idx - 1),
                    dp_matrix.get_insert(target_idx - 1, profile_idx - 1)
                        + profile.transition_score(Profile::INSERT_TO_MATCH_IDX, profile_idx - 1),
                    dp_matrix.get_special(target_idx - 1, Profile::SPECIAL_B_IDX)
                        + profile.transition_score(Profile::BEGIN_TO_MATCH_IDX, profile_idx - 1),
                    dp_matrix.get_delete(target_idx - 1, profile_idx - 1)
                        + profile.transition_score(Profile::DELETE_TO_MATCH_IDX, profile_idx - 1)
                ) + profile.match_score(current_target_character as usize, profile_idx),
            );

            // insert state
            dp_matrix.set_insert(
                target_idx,
                profile_idx,
                log_sum!(
                    dp_matrix.get_match(target_idx - 1, profile_idx)
                        + profile.transition_score(Profile::MATCH_TO_INSERT_IDX, profile_idx),
                    dp_matrix.get_insert(target_idx - 1, profile_idx)
                        + profile.transition_score(Profile::INSERT_TO_INSERT_IDX, profile_idx)
                ) + profile.insert_score(current_target_character as usize, profile_idx),
            );

            // delete state
            dp_matrix.set_delete(
                target_idx,
                profile_idx,
                log_sum!(
                    dp_matrix.get_match(target_idx, profile_idx - 1)
                        + profile.transition_score(Profile::MATCH_TO_DELETE_IDX, profile_idx - 1),
                    dp_matrix.get_delete(target_idx, profile_idx - 1)
                        + profile.transition_score(Profile::DELETE_TO_DELETE_IDX, profile_idx - 1)
                ),
            );

            // E state
            dp_matrix.set_special(
                target_idx,
                Profile::SPECIAL_E_IDX,
                log_sum!(
                    dp_matrix.get_match(target_idx, profile_idx) + end_score,
                    dp_matrix.get_delete(target_idx, profile_idx) + end_score,
                    dp_matrix.get_special(target_idx, Profile::SPECIAL_E_IDX)
                ),
            );
        }

        let last_profile_idx = params.right_row_bounds[target_idx];

        // unrolled match state match[M]
        dp_matrix.set_match(
            target_idx,
            last_profile_idx,
            log_sum!(
                dp_matrix.get_match(target_idx - 1, last_profile_idx - 1)
                    + profile.transition_score(Profile::MATCH_TO_MATCH_IDX, last_profile_idx - 1),
                dp_matrix.get_insert(target_idx - 1, last_profile_idx - 1)
                    + profile.transition_score(Profile::INSERT_TO_MATCH_IDX, last_profile_idx - 1),
                dp_matrix.get_special(target_idx - 1, Profile::SPECIAL_B_IDX)
                    + profile.transition_score(Profile::BEGIN_TO_MATCH_IDX, last_profile_idx - 1),
                dp_matrix.get_delete(target_idx - 1, last_profile_idx - 1)
                    + profile.transition_score(Profile::DELETE_TO_MATCH_IDX, last_profile_idx - 1)
            ) + profile.match_score(current_target_character as usize, last_profile_idx),
        );

        // unrolled insert state insert[M]
        dp_matrix.set_insert(target_idx, last_profile_idx, -f32::INFINITY);

        // unrolled delete state delete[M]
        dp_matrix.set_delete(
            target_idx,
            last_profile_idx,
            log_sum!(
                dp_matrix.get_match(target_idx, last_profile_idx - 1)
                    + profile.transition_score(Profile::MATCH_TO_DELETE_IDX, last_profile_idx - 1),
                dp_matrix.get_delete(target_idx, last_profile_idx - 1)
                    + profile.transition_score(Profile::DELETE_TO_DELETE_IDX, last_profile_idx - 1)
            ),
        );

        // unrolled E state
        dp_matrix.set_special(
            target_idx,
            Profile::SPECIAL_E_IDX,
            log_sum!(
                dp_matrix.get_match(target_idx, last_profile_idx),
                dp_matrix.get_delete(target_idx, last_profile_idx),
                dp_matrix.get_special(target_idx, Profile::SPECIAL_E_IDX)
            ),
        );

        // unrolled J state
        dp_matrix.set_special(
            target_idx,
            Profile::SPECIAL_J_IDX,
            log_sum!(
                dp_matrix.get_special(target_idx - 1, Profile::SPECIAL_J_IDX)
                    + profile.special_transition_score(
                        Profile::SPECIAL_J_IDX,
                        Profile::SPECIAL_LOOP_IDX
                    ),
                dp_matrix.get_special(target_idx, Profile::SPECIAL_E_IDX)
                    + profile.special_transition_score(
                        Profile::SPECIAL_E_IDX,
                        Profile::SPECIAL_LOOP_IDX
                    )
            ),
        );

        // unrolled C state
        dp_matrix.set_special(
            target_idx,
            Profile::SPECIAL_C_IDX,
            log_sum!(
                dp_matrix.get_special(target_idx - 1, Profile::SPECIAL_C_IDX)
                    + profile.special_transition_score(
                        Profile::SPECIAL_C_IDX,
                        Profile::SPECIAL_LOOP_IDX
                    ),
                dp_matrix.get_special(target_idx, Profile::SPECIAL_E_IDX)
                    + profile.special_transition_score(
                        Profile::SPECIAL_E_IDX,
                        Profile::SPECIAL_MOVE_IDX
                    )
            ),
        );

        // unrolled N state
        dp_matrix.set_special(
            target_idx,
            Profile::SPECIAL_N_IDX,
            dp_matrix.get_special(target_idx - 1, Profile::SPECIAL_N_IDX)
                + profile
                    .special_transition_score(Profile::SPECIAL_N_IDX, Profile::SPECIAL_LOOP_IDX),
        );

        // unrolled B state
        dp_matrix.set_special(
            target_idx,
            Profile::SPECIAL_B_IDX,
            log_sum!(
                dp_matrix.get_special(target_idx, Profile::SPECIAL_N_IDX)
                    + profile.special_transition_score(
                        Profile::SPECIAL_N_IDX,
                        Profile::SPECIAL_MOVE_IDX
                    ),
                dp_matrix.get_special(target_idx, Profile::SPECIAL_J_IDX)
                    + profile.special_transition_score(
                        Profile::SPECIAL_J_IDX,
                        Profile::SPECIAL_MOVE_IDX
                    )
            ),
        );
    }
}
