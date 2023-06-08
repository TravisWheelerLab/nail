use crate::log_sum;
use crate::structs::dp_matrix::DpMatrix;
use crate::structs::{Profile, Sequence};
use crate::util::log_add;

pub fn backward(profile: &Profile, target: &Sequence, dp_matrix: &mut impl DpMatrix) {
    let end_score: f32 = 0.0;

    // initialize the L row
    dp_matrix.set_special(target.length, Profile::SPECIAL_J_IDX, -f32::INFINITY);
    dp_matrix.set_special(target.length, Profile::SPECIAL_B_IDX, -f32::INFINITY);
    dp_matrix.set_special(target.length, Profile::SPECIAL_N_IDX, -f32::INFINITY);
    dp_matrix.set_special(
        target.length,
        Profile::SPECIAL_C_IDX,
        profile.special_transition_score(Profile::SPECIAL_C_IDX, Profile::SPECIAL_MOVE_IDX),
    );

    dp_matrix.set_special(
        target.length,
        Profile::SPECIAL_E_IDX,
        profile.special_transition_score(Profile::SPECIAL_C_IDX, Profile::SPECIAL_MOVE_IDX)
            + profile.special_transition_score(Profile::SPECIAL_E_IDX, Profile::SPECIAL_MOVE_IDX),
    );

    dp_matrix.set_match(
        target.length,
        profile.length,
        dp_matrix.get_special(target.length, Profile::SPECIAL_E_IDX),
    );

    dp_matrix.set_delete(
        target.length,
        profile.length,
        dp_matrix.get_special(target.length, Profile::SPECIAL_E_IDX),
    );

    dp_matrix.set_insert(target.length, profile.length, -f32::INFINITY);

    for profile_idx in (1..profile.length).rev() {
        dp_matrix.set_match(
            target.length,
            profile_idx,
            log_sum!(
                dp_matrix.get_special(target.length, Profile::SPECIAL_E_IDX) + end_score,
                dp_matrix.get_delete(target.length, profile_idx + 1)
                    + profile.transition_score(Profile::MATCH_TO_DELETE_IDX, profile_idx)
            ),
        );

        dp_matrix.set_insert(target.length, profile_idx, -f32::INFINITY);
        dp_matrix.set_delete(
            target.length,
            profile_idx,
            log_sum!(
                dp_matrix.get_special(target.length, Profile::SPECIAL_E_IDX) + end_score,
                dp_matrix.get_delete(target.length, profile_idx + 1)
                    + profile.transition_score(Profile::DELETE_TO_DELETE_IDX, profile_idx)
            ),
        );
    }

    // main recursion
    for target_idx in (1..target.length).rev() {
        let current_residue = target.digital_bytes[target_idx + 1] as usize;
        dp_matrix.set_special(
            target_idx,
            Profile::SPECIAL_B_IDX,
            dp_matrix.get_match(target_idx + 1, 1)
                + profile.transition_score(Profile::BEGIN_TO_MATCH_IDX, 0)
                + profile.match_score(current_residue, 1),
        );

        // note, this loops over the length of the profile, but it
        // incrementally changes the B state at index i (current target residue)
        for profile_idx in 2..=profile.length {
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
            profile.length,
            dp_matrix.get_special(target_idx, Profile::SPECIAL_E_IDX),
        );

        dp_matrix.set_insert(target_idx, profile.length, -f32::INFINITY);

        dp_matrix.set_delete(
            target_idx,
            profile.length,
            dp_matrix.get_special(target_idx, Profile::SPECIAL_E_IDX),
        );

        for profile_idx in (1..profile.length).rev() {
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

    let first_residue = target.digital_bytes[1] as usize;

    dp_matrix.set_special(
        0,
        Profile::SPECIAL_B_IDX,
        dp_matrix.get_match(1, 1)
            + profile.transition_score(Profile::BEGIN_TO_MATCH_IDX, 0)
            + profile.match_score(first_residue, 1),
    );

    for profile_idx in 2..=profile.length {
        dp_matrix.set_special(
            0,
            Profile::SPECIAL_B_IDX,
            log_sum!(
                dp_matrix.get_special(0, Profile::SPECIAL_B_IDX),
                dp_matrix.get_match(1, profile_idx)
                    + profile.transition_score(Profile::BEGIN_TO_MATCH_IDX, profile_idx - 1)
                    + profile.match_score(first_residue, profile_idx)
            ),
        );
    }

    dp_matrix.set_special(0, Profile::SPECIAL_J_IDX, -f32::INFINITY);
    dp_matrix.set_special(0, Profile::SPECIAL_C_IDX, -f32::INFINITY);
    dp_matrix.set_special(0, Profile::SPECIAL_E_IDX, -f32::INFINITY);
    dp_matrix.set_special(
        0,
        Profile::SPECIAL_N_IDX,
        log_sum!(
            dp_matrix.get_special(1, Profile::SPECIAL_N_IDX)
                + profile
                    .special_transition_score(Profile::SPECIAL_N_IDX, Profile::SPECIAL_LOOP_IDX),
            dp_matrix.get_special(0, Profile::SPECIAL_B_IDX)
                + profile
                    .special_transition_score(Profile::SPECIAL_N_IDX, Profile::SPECIAL_MOVE_IDX)
        ),
    );
    for profile_idx in (1..=profile.length).rev() {
        dp_matrix.set_match(0, profile_idx, -f32::INFINITY);
        dp_matrix.set_insert(0, profile_idx, -f32::INFINITY);
        dp_matrix.set_delete(0, profile_idx, -f32::INFINITY);
    }
}
