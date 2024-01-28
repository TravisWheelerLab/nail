use crate::align::structs::{DpMatrix, Trace};
use crate::structs::Profile;

pub fn traceback(
    profile: &Profile,
    posterior_matrix: &impl DpMatrix,
    optimal_matrix: &impl DpMatrix,
    trace: &mut Trace,
    target_end: usize,
) {
    let mut target_idx = target_end;
    let mut profile_idx = 0;

    let mut current_state_posterior_probability: f32;
    // we trace back starting from the last C state
    let mut previous_state: usize = Trace::C_STATE;
    let mut current_state: usize;

    trace.append_with_posterior_probability(Trace::T_STATE, target_idx, profile_idx, 0.0);
    trace.append_with_posterior_probability(Trace::C_STATE, target_idx, profile_idx, 0.0);

    while previous_state != Trace::S_STATE {
        current_state = match previous_state {
            Trace::C_STATE => {
                let c_to_c_path = profile.special_transition_score_delta(
                    Profile::SPECIAL_C_IDX,
                    Profile::SPECIAL_LOOP_IDX,
                ) * (optimal_matrix
                    .get_special(target_idx - 1, Profile::SPECIAL_C_IDX)
                    // TODO: why does this specific path involve a posterior probability?
                    + posterior_matrix.get_special(target_idx, Profile::SPECIAL_C_IDX));

                let c_to_e_path = profile
                    .transition_score_delta(Profile::SPECIAL_E_IDX, Profile::SPECIAL_MOVE_IDX)
                    * optimal_matrix.get_special(target_idx, Profile::SPECIAL_E_IDX);

                if c_to_c_path > c_to_e_path {
                    Trace::C_STATE
                } else {
                    Trace::E_STATE
                }
            }
            Trace::E_STATE => {
                let mut max_score = -f32::INFINITY;
                let mut state_of_max_score = 0;
                let mut profile_idx_of_max_score = 0;

                // TODO: why do we do this instead of just taking
                //       the state that follows the E state?
                for profile_idx in 1..=profile.length {
                    if optimal_matrix.get_match(target_idx, profile_idx) >= max_score {
                        max_score = optimal_matrix.get_match(target_idx, profile_idx);
                        state_of_max_score = Trace::M_STATE;
                        profile_idx_of_max_score = profile_idx;
                    }
                    if optimal_matrix.get_delete(target_idx, profile_idx) > max_score {
                        max_score = optimal_matrix.get_delete(target_idx, profile_idx);
                        state_of_max_score = Trace::D_STATE;
                        profile_idx_of_max_score = profile_idx;
                    }
                }
                profile_idx = profile_idx_of_max_score;
                state_of_max_score
            }
            Trace::M_STATE => {
                let possible_states: [usize; 4] = [
                    Trace::M_STATE,
                    Trace::I_STATE,
                    Trace::D_STATE,
                    Trace::B_STATE,
                ];

                let possible_paths: [f32; 4] = [
                    profile.transition_score_delta(Profile::MATCH_TO_MATCH_IDX, profile_idx - 1)
                        * optimal_matrix.get_match(target_idx - 1, profile_idx - 1),
                    profile.transition_score_delta(Profile::INSERT_TO_MATCH_IDX, profile_idx - 1)
                        * optimal_matrix.get_insert(target_idx - 1, profile_idx - 1),
                    profile.transition_score_delta(Profile::DELETE_TO_MATCH_IDX, profile_idx - 1)
                        * optimal_matrix.get_delete(target_idx - 1, profile_idx - 1),
                    profile.transition_score_delta(Profile::BEGIN_TO_MATCH_IDX, profile_idx - 1)
                        * optimal_matrix.get_special(target_idx - 1, Profile::SPECIAL_B_IDX),
                ];

                let mut argmax: usize = 0;
                for i in 1..4 {
                    if possible_paths[i] > possible_paths[argmax] {
                        argmax = i;
                    }
                }

                // a match means we have moved forward in the both the profile and the target
                profile_idx -= 1;
                target_idx -= 1;

                possible_states[argmax]
            }
            Trace::I_STATE => {
                let match_to_insert_path = profile
                    .transition_score_delta(Profile::MATCH_TO_INSERT_IDX, profile_idx)
                    * optimal_matrix.get_match(target_idx - 1, profile_idx);

                let insert_to_insert_path: f32 = profile
                    .transition_score_delta(Profile::INSERT_TO_INSERT_IDX, profile_idx)
                    * optimal_matrix.get_insert(target_idx - 1, profile_idx);

                // an insert means we moved forward only in the profile
                target_idx -= 1;

                if match_to_insert_path >= insert_to_insert_path {
                    Trace::M_STATE
                } else {
                    Trace::I_STATE
                }
            }
            Trace::D_STATE => {
                let match_to_delete_path = profile
                    .transition_score_delta(Profile::MATCH_TO_DELETE_IDX, profile_idx - 1)
                    * optimal_matrix.get_match(target_idx, profile_idx - 1);

                let delete_to_delete_path = profile
                    .transition_score_delta(Profile::DELETE_TO_DELETE_IDX, profile_idx - 1)
                    * optimal_matrix.get_delete(target_idx, profile_idx - 1);

                // a delete means we moved forward only in the profile
                profile_idx -= 1;

                if match_to_delete_path >= delete_to_delete_path {
                    Trace::M_STATE
                } else {
                    Trace::D_STATE
                }
            }
            Trace::B_STATE => {
                let n_to_b_path = profile.special_transition_score_delta(
                    Profile::SPECIAL_N_IDX,
                    Profile::SPECIAL_MOVE_IDX,
                ) * optimal_matrix
                    .get_special(target_idx, Profile::SPECIAL_N_IDX);

                let j_to_b_path = profile.special_transition_score_delta(
                    Profile::SPECIAL_J_IDX,
                    Profile::SPECIAL_MOVE_IDX,
                ) * optimal_matrix
                    .get_special(target_idx, Profile::SPECIAL_J_IDX);

                if n_to_b_path >= j_to_b_path {
                    Trace::N_STATE
                } else {
                    Trace::J_STATE
                }
            }
            Trace::N_STATE => {
                if target_idx == 0 {
                    Trace::S_STATE
                } else {
                    Trace::N_STATE
                }
            }
            Trace::J_STATE => {
                let j_to_j_path = profile.special_transition_score_delta(
                    Profile::SPECIAL_J_IDX,
                    Profile::SPECIAL_LOOP_IDX,
                ) * (optimal_matrix
                    .get_special(target_idx - 1, Profile::SPECIAL_J_IDX)
                    // TODO: why does this specific path involve a posterior probability?
                    + posterior_matrix.get_special(target_idx, Profile::SPECIAL_J_IDX));

                let e_to_j_path = profile.special_transition_score_delta(
                    Profile::SPECIAL_E_IDX,
                    Profile::SPECIAL_LOOP_IDX,
                ) * optimal_matrix
                    .get_special(target_idx, Profile::SPECIAL_E_IDX);

                if j_to_j_path > e_to_j_path {
                    Trace::J_STATE
                } else {
                    Trace::E_STATE
                }
            }
            _ => {
                panic!("bad state in traceback")
            }
        };

        current_state_posterior_probability = get_posterior_probability(
            posterior_matrix,
            current_state,
            previous_state,
            profile_idx,
            target_idx,
        );

        trace.append_with_posterior_probability(
            current_state,
            target_idx,
            profile_idx,
            current_state_posterior_probability,
        );

        if (current_state == Trace::N_STATE
            || current_state == Trace::J_STATE
            || current_state == Trace::C_STATE)
            && current_state == previous_state
        {
            target_idx -= 1;
        }
        previous_state = current_state;
    }
    trace.reverse();
}

pub fn get_posterior_probability(
    optimal_matrix: &impl DpMatrix,
    current_state: usize,
    previous_state: usize,
    profile_idx: usize,
    target_idx: usize,
) -> f32 {
    match current_state {
        Trace::M_STATE => optimal_matrix.get_match(target_idx, profile_idx),
        Trace::I_STATE => optimal_matrix.get_insert(target_idx, profile_idx),
        Trace::N_STATE => {
            if current_state == previous_state {
                optimal_matrix.get_special(target_idx, Profile::SPECIAL_N_IDX)
            } else {
                0.0
            }
        }
        Trace::C_STATE => {
            if current_state == previous_state {
                optimal_matrix.get_special(target_idx, Profile::SPECIAL_C_IDX)
            } else {
                0.0
            }
        }
        Trace::J_STATE => {
            if current_state == previous_state {
                optimal_matrix.get_special(target_idx, Profile::SPECIAL_J_IDX)
            } else {
                0.0
            }
        }
        _ => 0.0,
    }
}
