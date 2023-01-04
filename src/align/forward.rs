use crate::structs::profile::constants::{
    PROFILE_BEGIN_TO_MATCH, PROFILE_DELETE_TO_DELETE, PROFILE_DELETE_TO_MATCH,
    PROFILE_INSERT_TO_INSERT, PROFILE_INSERT_TO_MATCH, PROFILE_MATCH_TO_DELETE,
    PROFILE_MATCH_TO_INSERT, PROFILE_MATCH_TO_MATCH, SPECIAL_B, SPECIAL_C, SPECIAL_E, SPECIAL_J,
    SPECIAL_LOOP, SPECIAL_MOVE, SPECIAL_N,
};
use crate::structs::{DpMatrix, Profile, Sequence};
use crate::util::log_sum;
pub fn forward(profile: &Profile, target: &Sequence, dp_matrix: &mut DpMatrix) {
    let esc: f32 = 0.0;

    // initialize the zero row
    dp_matrix.set_special(0, SPECIAL_N, 0.0);
    dp_matrix.set_special(
        0,
        SPECIAL_B,
        profile.special_transition_score(SPECIAL_N, SPECIAL_MOVE),
    );
    dp_matrix.set_special(0, SPECIAL_E, -f32::INFINITY);
    dp_matrix.set_special(0, SPECIAL_C, -f32::INFINITY);
    dp_matrix.set_special(0, SPECIAL_J, -f32::INFINITY);

    for profile_idx in 0..=profile.length {
        dp_matrix.set_insert(0, profile_idx, -f32::INFINITY);
        dp_matrix.set_match(0, profile_idx, -f32::INFINITY);
        dp_matrix.set_delete(0, profile_idx, -f32::INFINITY);
    }

    // main recursion
    for target_idx in 1..=target.length {
        let current_residue = target.data[target_idx];
        let mut tmp_score: f32;

        dp_matrix.set_insert(target_idx, 0, -f32::INFINITY);
        dp_matrix.set_match(target_idx, 0, -f32::INFINITY);
        dp_matrix.set_delete(target_idx, 0, -f32::INFINITY);
        dp_matrix.set_special(target_idx, SPECIAL_E, -f32::INFINITY);

        for profile_idx in 1..profile.length {
            // match state
            tmp_score = log_sum(
                log_sum(
                    dp_matrix.get_match(target_idx - 1, profile_idx - 1)
                        + profile.transition_score(PROFILE_MATCH_TO_MATCH, profile_idx - 1),
                    dp_matrix.get_insert(target_idx - 1, profile_idx - 1)
                        + profile.transition_score(PROFILE_INSERT_TO_MATCH, profile_idx - 1),
                ),
                log_sum(
                    dp_matrix.get_special(target_idx - 1, SPECIAL_B)
                        + profile.transition_score(PROFILE_BEGIN_TO_MATCH, profile_idx - 1),
                    dp_matrix.get_delete(target_idx - 1, profile_idx - 1)
                        + profile.transition_score(PROFILE_DELETE_TO_MATCH, profile_idx - 1),
                ),
            );

            dp_matrix.set_match(
                target_idx,
                profile_idx,
                tmp_score + profile.match_score(current_residue as usize, profile_idx),
            );

            tmp_score = log_sum(
                dp_matrix.get_match(target_idx - 1, profile_idx)
                    + profile.transition_score(PROFILE_MATCH_TO_INSERT, profile_idx),
                dp_matrix.get_insert(target_idx - 1, profile_idx)
                    + profile.transition_score(PROFILE_INSERT_TO_INSERT, profile_idx),
            );

            // insert state
            dp_matrix.set_insert(
                target_idx,
                profile_idx,
                tmp_score + profile.insert_score(current_residue as usize, profile_idx),
            );

            // delete state
            dp_matrix.set_delete(
                target_idx,
                profile_idx,
                log_sum(
                    dp_matrix.get_match(target_idx, profile_idx - 1)
                        + profile.transition_score(PROFILE_MATCH_TO_DELETE, profile_idx - 1),
                    dp_matrix.get_delete(target_idx, profile_idx - 1)
                        + profile.transition_score(PROFILE_DELETE_TO_DELETE, profile_idx - 1),
                ),
            );

            // E state
            dp_matrix.set_special(
                target_idx,
                SPECIAL_E,
                log_sum(
                    log_sum(
                        dp_matrix.get_match(target_idx, profile_idx) + esc,
                        dp_matrix.get_delete(target_idx, profile_idx) + esc,
                    ),
                    dp_matrix.get_special(target_idx, SPECIAL_E),
                ),
            );
        }

        // unrolled match state match[M]
        tmp_score = log_sum(
            log_sum(
                dp_matrix.get_match(target_idx - 1, profile.length - 1)
                    + profile.transition_score(PROFILE_MATCH_TO_MATCH, profile.length - 1),
                dp_matrix.get_insert(target_idx - 1, profile.length - 1)
                    + profile.transition_score(PROFILE_INSERT_TO_MATCH, profile.length - 1),
            ),
            log_sum(
                dp_matrix.get_special(target_idx - 1, SPECIAL_B)
                    + profile.transition_score(PROFILE_BEGIN_TO_MATCH, profile.length - 1),
                dp_matrix.get_delete(target_idx - 1, profile.length - 1)
                    + profile.transition_score(PROFILE_DELETE_TO_MATCH, profile.length - 1),
            ),
        );

        dp_matrix.set_match(
            target_idx,
            profile.length,
            tmp_score + profile.match_score(current_residue as usize, profile.length),
        );

        // unrolled insert state insert[M]
        dp_matrix.set_insert(target_idx, profile.length, -f32::INFINITY);

        // unrolled delete state delete[M]
        dp_matrix.set_delete(
            target_idx,
            profile.length,
            log_sum(
                dp_matrix.get_match(target_idx, profile.length - 1)
                    + profile.transition_score(PROFILE_MATCH_TO_DELETE, profile.length - 1),
                dp_matrix.get_delete(target_idx, profile.length - 1)
                    + profile.transition_score(PROFILE_DELETE_TO_DELETE, profile.length - 1),
            ),
        );

        // unrolled E state
        dp_matrix.set_special(
            target_idx,
            SPECIAL_E,
            log_sum(
                log_sum(
                    dp_matrix.get_match(target_idx, profile.length),
                    dp_matrix.get_delete(target_idx, profile.length),
                ),
                dp_matrix.get_special(target_idx, SPECIAL_E),
            ),
        );

        // unrolled J state
        dp_matrix.set_special(
            target_idx,
            SPECIAL_J,
            log_sum(
                dp_matrix.get_special(target_idx - 1, SPECIAL_J)
                    + profile.special_transition_score(SPECIAL_J, SPECIAL_LOOP),
                dp_matrix.get_special(target_idx, SPECIAL_E)
                    + profile.special_transition_score(SPECIAL_E, SPECIAL_LOOP),
            ),
        );

        // unrolled C state
        dp_matrix.set_special(
            target_idx,
            SPECIAL_C,
            log_sum(
                dp_matrix.get_special(target_idx - 1, SPECIAL_C)
                    + profile.special_transition_score(SPECIAL_C, SPECIAL_LOOP),
                dp_matrix.get_special(target_idx, SPECIAL_E)
                    + profile.special_transition_score(SPECIAL_E, SPECIAL_MOVE),
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
            log_sum(
                dp_matrix.get_special(target_idx, SPECIAL_N)
                    + profile.special_transition_score(SPECIAL_N, SPECIAL_MOVE),
                dp_matrix.get_special(target_idx, SPECIAL_J)
                    + profile.special_transition_score(SPECIAL_J, SPECIAL_MOVE),
            ),
        );
    }
}
