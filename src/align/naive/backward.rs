use crate::log_sum;
use crate::structs::dp_matrix::DpMatrix;
use crate::structs::{Profile, Sequence};
use crate::timing::time;
use crate::util::log_add;
use anyhow::Result;

#[funci::timed(timer = time)]
pub fn backward(profile: &Profile, target: &Sequence, dp_matrix: &mut impl DpMatrix) -> Result<()> {
    let end_score: f32 = 0.0;

    // initialize the L row
    dp_matrix.set_special(target.length, Profile::SPECIAL_J, -f32::INFINITY);
    dp_matrix.set_special(target.length, Profile::SPECIAL_B, -f32::INFINITY);
    dp_matrix.set_special(target.length, Profile::SPECIAL_N, -f32::INFINITY);
    dp_matrix.set_special(
        target.length,
        Profile::SPECIAL_C,
        profile.special_transition_score(Profile::SPECIAL_C, Profile::SPECIAL_MOVE),
    );

    dp_matrix.set_special(
        target.length,
        Profile::SPECIAL_E,
        profile.special_transition_score(Profile::SPECIAL_C, Profile::SPECIAL_MOVE)
            + profile.special_transition_score(Profile::SPECIAL_E, Profile::SPECIAL_MOVE),
    );

    dp_matrix.set_match(
        target.length,
        profile.length,
        dp_matrix.get_special(target.length, Profile::SPECIAL_E),
    );

    dp_matrix.set_delete(
        target.length,
        profile.length,
        dp_matrix.get_special(target.length, Profile::SPECIAL_E),
    );

    dp_matrix.set_insert(target.length, profile.length, -f32::INFINITY);

    for profile_idx in (1..profile.length).rev() {
        dp_matrix.set_match(
            target.length,
            profile_idx,
            log_sum!(
                dp_matrix.get_special(target.length, Profile::SPECIAL_E) + end_score,
                dp_matrix.get_delete(target.length, profile_idx + 1)
                    + profile.transition_score(Profile::PROFILE_MATCH_TO_DELETE, profile_idx)
            ),
        );

        dp_matrix.set_insert(target.length, profile_idx, -f32::INFINITY);
        dp_matrix.set_delete(
            target.length,
            profile_idx,
            log_sum!(
                dp_matrix.get_special(target.length, Profile::SPECIAL_E) + end_score,
                dp_matrix.get_delete(target.length, profile_idx + 1)
                    + profile.transition_score(Profile::PROFILE_DELETE_TO_DELETE, profile_idx)
            ),
        );
    }

    // main recursion
    for target_idx in (1..target.length).rev() {
        let current_residue = target.digital_bytes[target_idx + 1] as usize;
        dp_matrix.set_special(
            target_idx,
            Profile::SPECIAL_B,
            dp_matrix.get_match(target_idx + 1, 1)
                + profile.transition_score(Profile::PROFILE_BEGIN_TO_MATCH, 0)
                + profile.match_score(current_residue, 1),
        );

        // note, this loops over the length of the profile, but it
        // incrementally changes the B state at index i (current target residue)
        // TODO: what?
        for profile_idx in 2..=profile.length {
            dp_matrix.set_special(
                target_idx,
                Profile::SPECIAL_B,
                log_sum!(
                    dp_matrix.get_special(target_idx, Profile::SPECIAL_B),
                    dp_matrix.get_match(target_idx + 1, profile_idx)
                        + profile
                            .transition_score(Profile::PROFILE_BEGIN_TO_MATCH, profile_idx - 1)
                        + profile.match_score(current_residue, profile_idx)
                ),
            );
        }

        dp_matrix.set_special(
            target_idx,
            Profile::SPECIAL_J,
            log_sum!(
                dp_matrix.get_special(target_idx + 1, Profile::SPECIAL_J)
                    + profile.special_transition_score(Profile::SPECIAL_J, Profile::SPECIAL_LOOP),
                dp_matrix.get_special(target_idx, Profile::SPECIAL_B)
                    + profile.special_transition_score(Profile::SPECIAL_J, Profile::SPECIAL_MOVE)
            ),
        );

        dp_matrix.set_special(
            target_idx,
            Profile::SPECIAL_C,
            dp_matrix.get_special(target_idx + 1, Profile::SPECIAL_C)
                + profile.special_transition_score(Profile::SPECIAL_C, Profile::SPECIAL_LOOP),
        );

        dp_matrix.set_special(
            target_idx,
            Profile::SPECIAL_E,
            log_sum!(
                dp_matrix.get_special(target_idx, Profile::SPECIAL_J)
                    + profile.special_transition_score(Profile::SPECIAL_E, Profile::SPECIAL_LOOP),
                dp_matrix.get_special(target_idx, Profile::SPECIAL_C)
                    + profile.special_transition_score(Profile::SPECIAL_E, Profile::SPECIAL_MOVE)
            ),
        );

        dp_matrix.set_special(
            target_idx,
            Profile::SPECIAL_N,
            log_sum!(
                dp_matrix.get_special(target_idx + 1, Profile::SPECIAL_N)
                    + profile.special_transition_score(Profile::SPECIAL_N, Profile::SPECIAL_LOOP),
                dp_matrix.get_special(target_idx, Profile::SPECIAL_B)
                    + profile.special_transition_score(Profile::SPECIAL_N, Profile::SPECIAL_MOVE)
            ),
        );

        dp_matrix.set_match(
            target_idx,
            profile.length,
            dp_matrix.get_special(target_idx, Profile::SPECIAL_E),
        );
        dp_matrix.set_insert(target_idx, profile.length, -f32::INFINITY);
        dp_matrix.set_delete(
            target_idx,
            profile.length,
            dp_matrix.get_special(target_idx, Profile::SPECIAL_E),
        );

        for profile_idx in (1..profile.length).rev() {
            dp_matrix.set_match(
                target_idx,
                profile_idx,
                log_sum!(
                    dp_matrix.get_match(target_idx + 1, profile_idx + 1)
                        + profile.transition_score(Profile::PROFILE_MATCH_TO_MATCH, profile_idx)
                        + profile.match_score(current_residue, profile_idx + 1),
                    dp_matrix.get_insert(target_idx + 1, profile_idx)
                        + profile.transition_score(Profile::PROFILE_MATCH_TO_INSERT, profile_idx)
                        + profile.insert_score(current_residue, profile_idx),
                    dp_matrix.get_special(target_idx, Profile::SPECIAL_E) + end_score,
                    dp_matrix.get_delete(target_idx, profile_idx + 1)
                        + profile.transition_score(Profile::PROFILE_MATCH_TO_DELETE, profile_idx)
                ),
            );

            dp_matrix.set_insert(
                target_idx,
                profile_idx,
                log_sum!(
                    dp_matrix.get_match(target_idx + 1, profile_idx + 1)
                        + profile.transition_score(Profile::PROFILE_INSERT_TO_MATCH, profile_idx)
                        + profile.match_score(current_residue, profile_idx + 1),
                    dp_matrix.get_insert(target_idx + 1, profile_idx)
                        + profile.transition_score(Profile::PROFILE_INSERT_TO_INSERT, profile_idx)
                        + profile.insert_score(current_residue, profile_idx)
                ),
            );

            dp_matrix.set_delete(
                target_idx,
                profile_idx,
                log_sum!(
                    dp_matrix.get_match(target_idx + 1, profile_idx + 1)
                        + profile.transition_score(Profile::PROFILE_DELETE_TO_MATCH, profile_idx)
                        + profile.match_score(current_residue, profile_idx + 1),
                    dp_matrix.get_delete(target_idx, profile_idx + 1)
                        + profile.transition_score(Profile::PROFILE_DELETE_TO_DELETE, profile_idx),
                    dp_matrix.get_special(target_idx, Profile::SPECIAL_E) + end_score
                ),
            );
        }
    }

    let first_residue = target.digital_bytes[1] as usize;

    dp_matrix.set_special(
        0,
        Profile::SPECIAL_B,
        dp_matrix.get_match(1, 1)
            + profile.transition_score(Profile::PROFILE_BEGIN_TO_MATCH, 0)
            + profile.match_score(first_residue, 1),
    );

    for profile_idx in 2..=profile.length {
        dp_matrix.set_special(
            0,
            Profile::SPECIAL_B,
            log_sum!(
                dp_matrix.get_special(0, Profile::SPECIAL_B),
                dp_matrix.get_match(1, profile_idx)
                    + profile.transition_score(Profile::PROFILE_BEGIN_TO_MATCH, profile_idx - 1)
                    + profile.match_score(first_residue, profile_idx)
            ),
        );
    }

    dp_matrix.set_special(0, Profile::SPECIAL_J, -f32::INFINITY);
    dp_matrix.set_special(0, Profile::SPECIAL_C, -f32::INFINITY);
    dp_matrix.set_special(0, Profile::SPECIAL_E, -f32::INFINITY);
    dp_matrix.set_special(
        0,
        Profile::SPECIAL_N,
        log_sum!(
            dp_matrix.get_special(1, Profile::SPECIAL_N)
                + profile.special_transition_score(Profile::SPECIAL_N, Profile::SPECIAL_LOOP),
            dp_matrix.get_special(0, Profile::SPECIAL_B)
                + profile.special_transition_score(Profile::SPECIAL_N, Profile::SPECIAL_MOVE)
        ),
    );
    for profile_idx in (1..=profile.length).rev() {
        dp_matrix.set_match(0, profile_idx, -f32::INFINITY);
        dp_matrix.set_insert(0, profile_idx, -f32::INFINITY);
        dp_matrix.set_delete(0, profile_idx, -f32::INFINITY);
    }

    Ok(())
}
