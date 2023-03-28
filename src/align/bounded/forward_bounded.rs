use crate::align::bounded::structs::RowBoundParams;
use crate::log_sum;
use crate::structs::dp_matrix::DpMatrix;
use crate::structs::profile::constants::{
    PROFILE_BEGIN_TO_MATCH, PROFILE_DELETE_TO_DELETE, PROFILE_DELETE_TO_MATCH,
    PROFILE_INSERT_TO_INSERT, PROFILE_INSERT_TO_MATCH, PROFILE_MATCH_TO_DELETE,
    PROFILE_MATCH_TO_INSERT, PROFILE_MATCH_TO_MATCH, SPECIAL_B, SPECIAL_C, SPECIAL_E, SPECIAL_J,
    SPECIAL_LOOP, SPECIAL_MOVE, SPECIAL_N,
};
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

    dp_matrix.set_special(params.target_start - 1, SPECIAL_N, 0.0);
    dp_matrix.set_special(
        params.target_start - 1,
        SPECIAL_B,
        profile.special_transition_score(SPECIAL_N, SPECIAL_MOVE),
    );
    dp_matrix.set_special(params.target_start - 1, SPECIAL_E, -f32::INFINITY);
    dp_matrix.set_special(params.target_start - 1, SPECIAL_C, -f32::INFINITY);
    dp_matrix.set_special(params.target_start - 1, SPECIAL_J, -f32::INFINITY);

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
                        + profile.transition_score(PROFILE_MATCH_TO_MATCH, profile_idx - 1),
                    dp_matrix.get_insert(target_idx - 1, profile_idx - 1)
                        + profile.transition_score(PROFILE_INSERT_TO_MATCH, profile_idx - 1),
                    dp_matrix.get_special(target_idx - 1, SPECIAL_B)
                        + profile.transition_score(PROFILE_BEGIN_TO_MATCH, profile_idx - 1),
                    dp_matrix.get_delete(target_idx - 1, profile_idx - 1)
                        + profile.transition_score(PROFILE_DELETE_TO_MATCH, profile_idx - 1)
                ) + profile.match_score(current_target_character as usize, profile_idx),
            );

            // insert state
            dp_matrix.set_insert(
                target_idx,
                profile_idx,
                log_sum!(
                    dp_matrix.get_match(target_idx - 1, profile_idx)
                        + profile.transition_score(PROFILE_MATCH_TO_INSERT, profile_idx),
                    dp_matrix.get_insert(target_idx - 1, profile_idx)
                        + profile.transition_score(PROFILE_INSERT_TO_INSERT, profile_idx)
                ) + profile.insert_score(current_target_character as usize, profile_idx),
            );

            // delete state
            dp_matrix.set_delete(
                target_idx,
                profile_idx,
                log_sum!(
                    dp_matrix.get_match(target_idx, profile_idx - 1)
                        + profile.transition_score(PROFILE_MATCH_TO_DELETE, profile_idx - 1),
                    dp_matrix.get_delete(target_idx, profile_idx - 1)
                        + profile.transition_score(PROFILE_DELETE_TO_DELETE, profile_idx - 1)
                ),
            );

            // E state
            dp_matrix.set_special(
                target_idx,
                SPECIAL_E,
                log_sum!(
                    dp_matrix.get_match(target_idx, profile_idx) + end_score,
                    dp_matrix.get_delete(target_idx, profile_idx) + end_score,
                    dp_matrix.get_special(target_idx, SPECIAL_E)
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
                    + profile.transition_score(PROFILE_MATCH_TO_MATCH, last_profile_idx - 1),
                dp_matrix.get_insert(target_idx - 1, last_profile_idx - 1)
                    + profile.transition_score(PROFILE_INSERT_TO_MATCH, last_profile_idx - 1),
                dp_matrix.get_special(target_idx - 1, SPECIAL_B)
                    + profile.transition_score(PROFILE_BEGIN_TO_MATCH, last_profile_idx - 1),
                dp_matrix.get_delete(target_idx - 1, last_profile_idx - 1)
                    + profile.transition_score(PROFILE_DELETE_TO_MATCH, last_profile_idx - 1)
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
                    + profile.transition_score(PROFILE_MATCH_TO_DELETE, last_profile_idx - 1),
                dp_matrix.get_delete(target_idx, last_profile_idx - 1)
                    + profile.transition_score(PROFILE_DELETE_TO_DELETE, last_profile_idx - 1)
            ),
        );

        // unrolled E state
        dp_matrix.set_special(
            target_idx,
            SPECIAL_E,
            log_sum!(
                dp_matrix.get_match(target_idx, last_profile_idx),
                dp_matrix.get_delete(target_idx, last_profile_idx),
                dp_matrix.get_special(target_idx, SPECIAL_E)
            ),
        );

        // unrolled J state
        dp_matrix.set_special(
            target_idx,
            SPECIAL_J,
            log_sum!(
                dp_matrix.get_special(target_idx - 1, SPECIAL_J)
                    + profile.special_transition_score(SPECIAL_J, SPECIAL_LOOP),
                dp_matrix.get_special(target_idx, SPECIAL_E)
                    + profile.special_transition_score(SPECIAL_E, SPECIAL_LOOP)
            ),
        );

        // unrolled C state
        dp_matrix.set_special(
            target_idx,
            SPECIAL_C,
            log_sum!(
                dp_matrix.get_special(target_idx - 1, SPECIAL_C)
                    + profile.special_transition_score(SPECIAL_C, SPECIAL_LOOP),
                dp_matrix.get_special(target_idx, SPECIAL_E)
                    + profile.special_transition_score(SPECIAL_E, SPECIAL_MOVE)
            ),
        );

        // unrolled N state
        dp_matrix.set_special(
            target_idx,
            SPECIAL_N,
            dp_matrix.get_special(target_idx - 1, SPECIAL_N)
                + profile.special_transition_score(SPECIAL_N, SPECIAL_LOOP),
        );

        // unrolled B state
        dp_matrix.set_special(
            target_idx,
            SPECIAL_B,
            log_sum!(
                dp_matrix.get_special(target_idx, SPECIAL_N)
                    + profile.special_transition_score(SPECIAL_N, SPECIAL_MOVE),
                dp_matrix.get_special(target_idx, SPECIAL_J)
                    + profile.special_transition_score(SPECIAL_J, SPECIAL_MOVE)
            ),
        );
    }
}
