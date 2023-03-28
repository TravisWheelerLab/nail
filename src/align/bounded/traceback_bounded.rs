use crate::structs::dp_matrix::DpMatrix;
use crate::structs::profile::constants::{
    PROFILE_BEGIN_TO_MATCH, PROFILE_DELETE_TO_DELETE, PROFILE_DELETE_TO_MATCH,
    PROFILE_INSERT_TO_INSERT, PROFILE_INSERT_TO_MATCH, PROFILE_MATCH_TO_DELETE,
    PROFILE_MATCH_TO_INSERT, PROFILE_MATCH_TO_MATCH, SPECIAL_B, SPECIAL_C, SPECIAL_E, SPECIAL_J,
    SPECIAL_LOOP, SPECIAL_MOVE, SPECIAL_N,
};
use crate::structs::trace::constants::{
    TRACE_B, TRACE_C, TRACE_D, TRACE_E, TRACE_I, TRACE_J, TRACE_M, TRACE_N, TRACE_S, TRACE_T,
};
use crate::structs::{Profile, Trace};

use crate::timing::time;

#[funci::timed(timer = time)]
pub fn traceback_bounded(
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
    let mut previous_state: usize = TRACE_C;
    let mut current_state: usize;

    trace.append_with_posterior_probability(TRACE_T, target_idx, profile_idx, 0.0);
    trace.append_with_posterior_probability(TRACE_C, target_idx, profile_idx, 0.0);

    while previous_state != TRACE_S {
        current_state = match previous_state {
            TRACE_C => {
                let c_to_c_path = profile.special_transition_score_delta(SPECIAL_C, SPECIAL_LOOP)
                    * (optimal_matrix.get_special(target_idx - 1, SPECIAL_C)
                        + posterior_matrix.get_special(target_idx, SPECIAL_C));

                let c_to_e_path = profile.transition_score_delta(SPECIAL_E, SPECIAL_MOVE)
                    * optimal_matrix.get_special(target_idx, SPECIAL_E);

                if c_to_c_path > c_to_e_path {
                    TRACE_C
                } else {
                    TRACE_E
                }
            }
            TRACE_E => {
                let mut max_score = -f32::INFINITY;
                let mut state_of_max_score = 0;
                let mut profile_idx_of_max_score = 0;

                for profile_idx in 1..=profile.length {
                    if optimal_matrix.get_match(target_idx, profile_idx) >= max_score {
                        max_score = optimal_matrix.get_match(target_idx, profile_idx);
                        state_of_max_score = TRACE_M;
                        profile_idx_of_max_score = profile_idx;
                    }
                    if optimal_matrix.get_delete(target_idx, profile_idx) > max_score {
                        max_score = optimal_matrix.get_delete(target_idx, profile_idx);
                        state_of_max_score = TRACE_D;
                        profile_idx_of_max_score = profile_idx;
                    }
                }
                profile_idx = profile_idx_of_max_score;
                state_of_max_score
            }
            TRACE_M => {
                let possible_states: [usize; 4] = [TRACE_M, TRACE_I, TRACE_D, TRACE_B];
                let possible_paths: [f32; 4] = [
                    profile.transition_score_delta(PROFILE_MATCH_TO_MATCH, profile_idx - 1)
                        * optimal_matrix.get_match(target_idx - 1, profile_idx - 1),
                    profile.transition_score_delta(PROFILE_INSERT_TO_MATCH, profile_idx - 1)
                        * optimal_matrix.get_insert(target_idx - 1, profile_idx - 1),
                    profile.transition_score_delta(PROFILE_DELETE_TO_MATCH, profile_idx - 1)
                        * optimal_matrix.get_delete(target_idx - 1, profile_idx - 1),
                    profile.transition_score_delta(PROFILE_BEGIN_TO_MATCH, profile_idx - 1)
                        * optimal_matrix.get_special(target_idx - 1, SPECIAL_B),
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
            TRACE_I => {
                let match_to_insert_path = profile
                    .transition_score_delta(PROFILE_MATCH_TO_INSERT, profile_idx)
                    * optimal_matrix.get_match(target_idx - 1, profile_idx);
                let insert_to_insert_path: f32 = profile
                    .transition_score_delta(PROFILE_INSERT_TO_INSERT, profile_idx)
                    * optimal_matrix.get_insert(target_idx - 1, profile_idx);

                // an insert means we moved forward only in the profile
                target_idx -= 1;

                if match_to_insert_path >= insert_to_insert_path {
                    TRACE_M
                } else {
                    TRACE_I
                }
            }
            TRACE_D => {
                let match_to_delete_path = profile
                    .transition_score_delta(PROFILE_MATCH_TO_DELETE, profile_idx - 1)
                    * optimal_matrix.get_match(target_idx, profile_idx - 1);
                let delete_to_delete_path = profile
                    .transition_score_delta(PROFILE_DELETE_TO_DELETE, profile_idx - 1)
                    * optimal_matrix.get_delete(target_idx, profile_idx - 1);

                // a delete means we moved forward only in the profile
                profile_idx -= 1;

                if match_to_delete_path >= delete_to_delete_path {
                    TRACE_M
                } else {
                    TRACE_D
                }
            }
            TRACE_B => {
                let n_to_b_path = profile.special_transition_score_delta(SPECIAL_N, SPECIAL_MOVE)
                    * optimal_matrix.get_match(target_idx, SPECIAL_N);
                let j_to_b_path = profile.transition_score_delta(SPECIAL_J, SPECIAL_MOVE)
                    * optimal_matrix.get_special(target_idx, SPECIAL_J);

                if n_to_b_path >= j_to_b_path {
                    TRACE_N
                } else {
                    TRACE_J
                }
            }
            TRACE_N => {
                if target_idx == 0 {
                    TRACE_S
                } else {
                    TRACE_N
                }
            }
            TRACE_J => {
                let j_to_j_path = profile.special_transition_score_delta(SPECIAL_J, SPECIAL_LOOP)
                    * (optimal_matrix.get_special(target_idx - 1, SPECIAL_J)
                        + posterior_matrix.get_special(target_idx, SPECIAL_J));
                let e_to_j_path = profile.transition_score_delta(SPECIAL_E, SPECIAL_LOOP)
                    * optimal_matrix.get_special(target_idx, SPECIAL_E);

                if j_to_j_path > e_to_j_path {
                    TRACE_J
                } else {
                    TRACE_E
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

        if (current_state == TRACE_N || current_state == TRACE_J || current_state == TRACE_C)
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
        TRACE_M => optimal_matrix.get_match(target_idx, profile_idx),
        TRACE_I => optimal_matrix.get_insert(target_idx, profile_idx),
        TRACE_N => {
            if current_state == previous_state {
                optimal_matrix.get_special(target_idx, SPECIAL_N)
            } else {
                0.0
            }
        }
        TRACE_C => {
            if current_state == previous_state {
                optimal_matrix.get_special(target_idx, SPECIAL_C)
            } else {
                0.0
            }
        }
        TRACE_J => {
            if current_state == previous_state {
                optimal_matrix.get_special(target_idx, SPECIAL_J)
            } else {
                0.0
            }
        }
        _ => 0.0,
    }
}
