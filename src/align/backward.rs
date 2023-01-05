use crate::structs::profile::constants::{
    PROFILE_BEGIN_TO_MATCH, PROFILE_DELETE_TO_DELETE, PROFILE_DELETE_TO_MATCH,
    PROFILE_INSERT_TO_INSERT, PROFILE_INSERT_TO_MATCH, PROFILE_MATCH_TO_DELETE,
    PROFILE_MATCH_TO_INSERT, PROFILE_MATCH_TO_MATCH, SPECIAL_B, SPECIAL_C, SPECIAL_E, SPECIAL_J,
    SPECIAL_LOOP, SPECIAL_MOVE, SPECIAL_N,
};
use crate::structs::{DpMatrix, Profile, Sequence};
use crate::util::log_sum;

pub fn backward(profile: &Profile, target: &Sequence, dp_matrix: &mut DpMatrix) {
    let esc: f32 = 0.0;

    // initialize the L row
    dp_matrix.set_special(target.length, SPECIAL_J, -f32::INFINITY);
    dp_matrix.set_special(target.length, SPECIAL_B, -f32::INFINITY);
    dp_matrix.set_special(target.length, SPECIAL_N, -f32::INFINITY);
    dp_matrix.set_special(
        target.length,
        SPECIAL_C,
        profile.special_transition_score(SPECIAL_C, SPECIAL_MOVE),
    );

    dp_matrix.set_special(
        target.length,
        SPECIAL_E,
        profile.special_transition_score(SPECIAL_C, SPECIAL_MOVE)
            + profile.special_transition_score(SPECIAL_E, SPECIAL_MOVE),
    );

    dp_matrix.set_match(
        target.length,
        profile.length,
        dp_matrix.get_special(target.length, SPECIAL_E),
    );

    dp_matrix.set_delete(
        target.length,
        profile.length,
        dp_matrix.get_special(target.length, SPECIAL_E),
    );

    dp_matrix.set_insert(target.length, profile.length, -f32::INFINITY);

    for profile_idx in (1..profile.length).rev() {
        dp_matrix.set_match(
            target.length,
            profile_idx,
            log_sum(
                dp_matrix.get_special(target.length, SPECIAL_E) + esc,
                dp_matrix.get_delete(target.length, profile_idx + 1)
                    + profile.transition_score(PROFILE_MATCH_TO_DELETE, profile_idx),
            ),
        );

        dp_matrix.set_insert(target.length, profile_idx, -f32::INFINITY);
        dp_matrix.set_delete(
            target.length,
            profile_idx,
            log_sum(
                dp_matrix.get_special(target.length, SPECIAL_E) + esc,
                dp_matrix.get_delete(target.length, profile_idx + 1)
                    + profile.transition_score(PROFILE_DELETE_TO_DELETE, profile_idx),
            ),
        );
    }

    // main recursion
    for target_idx in (1..target.length).rev() {
        let current_residue = target.digital_bytes[target_idx + 1] as usize;
        dp_matrix.set_special(
            target_idx,
            SPECIAL_B,
            dp_matrix.get_match(target_idx + 1, 1)
                + profile.transition_score(PROFILE_BEGIN_TO_MATCH, 0)
                + profile.match_score(current_residue as usize, 1),
        );

        // note, this loops over the length of the profile, but it
        // incrementally changes the B state at index i (current target residue)
        // TODO: what?
        for profile_idx in 2..=profile.length {
            dp_matrix.set_special(
                target_idx,
                SPECIAL_B,
                log_sum(
                    dp_matrix.get_special(target_idx, SPECIAL_B),
                    dp_matrix.get_match(target_idx + 1, profile_idx)
                        + profile.transition_score(PROFILE_BEGIN_TO_MATCH, profile_idx - 1)
                        + profile.match_score(current_residue, profile_idx),
                ),
            );
        }

        dp_matrix.set_special(
            target_idx,
            SPECIAL_J,
            log_sum(
                dp_matrix.get_special(target_idx + 1, SPECIAL_J)
                    + profile.special_transition_score(SPECIAL_J, SPECIAL_LOOP),
                dp_matrix.get_special(target_idx, SPECIAL_B)
                    + profile.special_transition_score(SPECIAL_J, SPECIAL_MOVE),
            ),
        );

        dp_matrix.set_special(
            target_idx,
            SPECIAL_C,
            dp_matrix.get_special(target_idx + 1, SPECIAL_C)
                + profile.special_transition_score(SPECIAL_C, SPECIAL_LOOP),
        );

        dp_matrix.set_special(
            target_idx,
            SPECIAL_E,
            log_sum(
                dp_matrix.get_special(target_idx, SPECIAL_J)
                    + profile.special_transition_score(SPECIAL_E, SPECIAL_LOOP),
                dp_matrix.get_special(target_idx, SPECIAL_C)
                    + profile.special_transition_score(SPECIAL_E, SPECIAL_MOVE),
            ),
        );

        dp_matrix.set_special(
            target_idx,
            SPECIAL_N,
            log_sum(
                dp_matrix.get_special(target_idx + 1, SPECIAL_N)
                    + profile.special_transition_score(SPECIAL_N, SPECIAL_LOOP),
                dp_matrix.get_special(target_idx, SPECIAL_B)
                    + profile.special_transition_score(SPECIAL_N, SPECIAL_MOVE),
            ),
        );

        dp_matrix.set_match(
            target_idx,
            profile.length,
            dp_matrix.get_special(target_idx, SPECIAL_E),
        );
        dp_matrix.set_insert(target_idx, profile.length, -f32::INFINITY);
        dp_matrix.set_delete(
            target_idx,
            profile.length,
            dp_matrix.get_special(target_idx, SPECIAL_E),
        );

        for profile_idx in (1..profile.length).rev() {
            dp_matrix.set_match(
                target_idx,
                profile_idx,
                log_sum(
                    log_sum(
                        dp_matrix.get_match(target_idx + 1, profile_idx + 1)
                            + profile.transition_score(PROFILE_MATCH_TO_MATCH, profile_idx)
                            + profile.match_score(current_residue, profile_idx + 1),
                        dp_matrix.get_insert(target_idx + 1, profile_idx)
                            + profile.transition_score(PROFILE_MATCH_TO_INSERT, profile_idx)
                            + profile.insert_score(current_residue, profile_idx),
                    ),
                    log_sum(
                        dp_matrix.get_special(target_idx, SPECIAL_E) + esc,
                        dp_matrix.get_delete(target_idx, profile_idx + 1)
                            + profile.transition_score(PROFILE_MATCH_TO_DELETE, profile_idx),
                    ),
                ),
            );

            dp_matrix.set_insert(
                target_idx,
                profile_idx,
                log_sum(
                    dp_matrix.get_match(target_idx + 1, profile_idx + 1)
                        + profile.transition_score(PROFILE_INSERT_TO_MATCH, profile_idx)
                        + profile.match_score(current_residue, profile_idx + 1),
                    dp_matrix.get_insert(target_idx + 1, profile_idx)
                        + profile.transition_score(PROFILE_INSERT_TO_INSERT, profile_idx)
                        + profile.insert_score(current_residue, profile_idx),
                ),
            );

            dp_matrix.set_delete(
                target_idx,
                profile_idx,
                log_sum(
                    dp_matrix.get_match(target_idx + 1, profile_idx + 1)
                        + profile.transition_score(PROFILE_DELETE_TO_MATCH, profile_idx)
                        + profile.match_score(current_residue, profile_idx + 1),
                    log_sum(
                        dp_matrix.get_delete(target_idx, profile_idx + 1)
                            + profile.transition_score(PROFILE_DELETE_TO_DELETE, profile_idx),
                        dp_matrix.get_special(target_idx, SPECIAL_E) + esc,
                    ),
                ),
            );
        }
    }

    let first_residue = target.digital_bytes[1] as usize;

    dp_matrix.set_special(
        0,
        SPECIAL_B,
        dp_matrix.get_match(1, 1)
            + profile.transition_score(PROFILE_BEGIN_TO_MATCH, 0)
            + profile.match_score(first_residue, 1),
    );

    for profile_idx in 2..=profile.length {
        dp_matrix.set_special(
            0,
            SPECIAL_B,
            log_sum(
                dp_matrix.get_special(0, SPECIAL_B),
                dp_matrix.get_match(1, profile_idx)
                    + profile.transition_score(PROFILE_BEGIN_TO_MATCH, profile_idx - 1)
                    + profile.match_score(first_residue, profile_idx),
            ),
        );
    }

    dp_matrix.set_special(0, SPECIAL_J, -f32::INFINITY);
    dp_matrix.set_special(0, SPECIAL_C, -f32::INFINITY);
    dp_matrix.set_special(0, SPECIAL_E, -f32::INFINITY);
    dp_matrix.set_special(
        0,
        SPECIAL_N,
        log_sum(
            dp_matrix.get_special(1, SPECIAL_N)
                + profile.special_transition_score(SPECIAL_N, SPECIAL_LOOP),
            dp_matrix.get_special(0, SPECIAL_B)
                + profile.special_transition_score(SPECIAL_N, SPECIAL_MOVE),
        ),
    );
    for profile_idx in (1..=profile.length).rev() {
        dp_matrix.set_match(0, profile_idx, -f32::INFINITY);
        dp_matrix.set_insert(0, profile_idx, -f32::INFINITY);
        dp_matrix.set_delete(0, profile_idx, -f32::INFINITY);
    }
}
