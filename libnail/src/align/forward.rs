use crate::align::structs::{DpMatrix, RowBounds};
use crate::log_sum;
use crate::structs::{Profile, Sequence};
use crate::util::log_add;

use super::Nats;

pub fn forward(
    profile: &Profile,
    target: &Sequence,
    dp_matrix: &mut impl DpMatrix,
    bounds: &RowBounds,
) -> Nats {
    let end_score: f32 = 0.0;

    dp_matrix.set_special(bounds.seq_start - 1, Profile::N_IDX, 0.0);

    dp_matrix.set_special(
        bounds.seq_start - 1,
        Profile::B_IDX,
        profile.special_transition_score(Profile::N_IDX, Profile::SPECIAL_MOVE_IDX),
    );
    dp_matrix.set_special(bounds.seq_start - 1, Profile::E_IDX, -f32::INFINITY);
    dp_matrix.set_special(bounds.seq_start - 1, Profile::C_IDX, -f32::INFINITY);

    for target_idx in bounds.seq_start..=bounds.seq_end {
        let current_target_character = target.digital_bytes[target_idx];

        for profile_idx in bounds.left_row_bounds[target_idx]..bounds.right_row_bounds[target_idx] {
            // match state
            dp_matrix.set_match(
                target_idx,
                profile_idx,
                log_sum!(
                    dp_matrix.get_match(target_idx - 1, profile_idx - 1)
                        + profile.transition_score(Profile::M_M_IDX, profile_idx - 1),
                    dp_matrix.get_insert(target_idx - 1, profile_idx - 1)
                        + profile.transition_score(Profile::I_M_IDX, profile_idx - 1),
                    dp_matrix.get_special(target_idx - 1, Profile::B_IDX)
                        + profile.transition_score(Profile::B_M_IDX, profile_idx - 1),
                    dp_matrix.get_delete(target_idx - 1, profile_idx - 1)
                        + profile.transition_score(Profile::D_M_IDX, profile_idx - 1)
                ) + profile.match_score(current_target_character as usize, profile_idx),
            );

            // insert state
            dp_matrix.set_insert(
                target_idx,
                profile_idx,
                log_sum!(
                    dp_matrix.get_match(target_idx - 1, profile_idx)
                        + profile.transition_score(Profile::M_I_IDX, profile_idx),
                    dp_matrix.get_insert(target_idx - 1, profile_idx)
                        + profile.transition_score(Profile::I_I_IDX, profile_idx)
                ) + profile.insert_score(current_target_character as usize, profile_idx),
            );

            // delete state
            dp_matrix.set_delete(
                target_idx,
                profile_idx,
                log_sum!(
                    dp_matrix.get_match(target_idx, profile_idx - 1)
                        + profile.transition_score(Profile::M_D_IDX, profile_idx - 1),
                    dp_matrix.get_delete(target_idx, profile_idx - 1)
                        + profile.transition_score(Profile::D_D_IDX, profile_idx - 1)
                ),
            );

            // E state
            dp_matrix.set_special(
                target_idx,
                Profile::E_IDX,
                log_sum!(
                    dp_matrix.get_match(target_idx, profile_idx) + end_score,
                    dp_matrix.get_delete(target_idx, profile_idx) + end_score,
                    dp_matrix.get_special(target_idx, Profile::E_IDX)
                ),
            );
        }

        let last_profile_idx = bounds.right_row_bounds[target_idx];

        // unrolled match state match[M]
        dp_matrix.set_match(
            target_idx,
            last_profile_idx,
            log_sum!(
                dp_matrix.get_match(target_idx - 1, last_profile_idx - 1)
                    + profile.transition_score(Profile::M_M_IDX, last_profile_idx - 1),
                dp_matrix.get_insert(target_idx - 1, last_profile_idx - 1)
                    + profile.transition_score(Profile::I_M_IDX, last_profile_idx - 1),
                dp_matrix.get_special(target_idx - 1, Profile::B_IDX)
                    + profile.transition_score(Profile::B_M_IDX, last_profile_idx - 1),
                dp_matrix.get_delete(target_idx - 1, last_profile_idx - 1)
                    + profile.transition_score(Profile::D_M_IDX, last_profile_idx - 1)
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
                    + profile.transition_score(Profile::M_D_IDX, last_profile_idx - 1),
                dp_matrix.get_delete(target_idx, last_profile_idx - 1)
                    + profile.transition_score(Profile::D_D_IDX, last_profile_idx - 1)
            ),
        );

        // unrolled E state
        dp_matrix.set_special(
            target_idx,
            Profile::E_IDX,
            log_sum!(
                dp_matrix.get_match(target_idx, last_profile_idx),
                dp_matrix.get_delete(target_idx, last_profile_idx),
                dp_matrix.get_special(target_idx, Profile::E_IDX)
            ),
        );

        // unrolled C state
        dp_matrix.set_special(
            target_idx,
            Profile::C_IDX,
            log_sum!(
                dp_matrix.get_special(target_idx - 1, Profile::C_IDX)
                    + profile.special_transition_score(Profile::C_IDX, Profile::SPECIAL_LOOP_IDX),
                dp_matrix.get_special(target_idx, Profile::E_IDX)
                    + profile.special_transition_score(Profile::E_IDX, Profile::SPECIAL_MOVE_IDX)
            ),
        );

        // unrolled N state
        dp_matrix.set_special(
            target_idx,
            Profile::N_IDX,
            dp_matrix.get_special(target_idx - 1, Profile::N_IDX)
                + profile.special_transition_score(Profile::N_IDX, Profile::SPECIAL_LOOP_IDX),
        );

        // unrolled B state
        dp_matrix.set_special(
            target_idx,
            Profile::B_IDX,
            log_sum!(
                dp_matrix.get_special(target_idx, Profile::N_IDX)
                    + profile.special_transition_score(Profile::N_IDX, Profile::SPECIAL_MOVE_IDX)
            ),
        );
    }

    // what: sum up the loop transitions to the N and/or C states
    //       once for every position in the target sequence that
    //       isn't included in the cloud
    // why: we want the best approximation of the full Forward score
    let aligned_target_length = bounds.seq_end - bounds.seq_start + 1;
    let unaligned_target_length = target.length - aligned_target_length;
    let background_correction = unaligned_target_length as f32
        * profile.special_transition_score(Profile::N_IDX, Profile::SPECIAL_LOOP_IDX);

    let final_c_state_score = dp_matrix.get_special(bounds.seq_end, Profile::C_IDX);
    let c_to_exit_score =
        profile.special_transition_score(Profile::C_IDX, Profile::SPECIAL_MOVE_IDX);

    Nats(final_c_state_score + background_correction + c_to_exit_score)
}
