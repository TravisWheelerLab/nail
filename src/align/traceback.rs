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
use crate::structs::{DpMatrix3D, Profile, Trace};

pub fn traceback(
    profile: &Profile,
    posterior_matrix: &DpMatrix3D,
    optimal_matrix: &DpMatrix3D,
    trace: &mut Trace,
) {
    let mut target_idx = optimal_matrix.target_length;
    let mut profile_idx = 0;

    let mut posterior_probability: f32;
    let mut previous_state: usize = TRACE_C;
    let mut current_state: usize;

    trace.append_with_posterior_probability(TRACE_T, target_idx, profile_idx, 0.0);
    trace.append_with_posterior_probability(TRACE_C, target_idx, profile_idx, 0.0);

    while previous_state != TRACE_S {
        match previous_state {
            TRACE_M => {
                current_state = select_m(profile, optimal_matrix, target_idx, profile_idx);
                profile_idx -= 1;
                target_idx -= 1;
            }
            TRACE_D => {
                current_state = select_d(profile, optimal_matrix, target_idx, profile_idx);
                profile_idx -= 1;
            }
            TRACE_I => {
                current_state = select_i(profile, optimal_matrix, target_idx, profile_idx);
                target_idx -= 1;
            }
            TRACE_N => {
                current_state = select_n(target_idx);
            }
            TRACE_C => {
                current_state = select_c(profile, posterior_matrix, optimal_matrix, target_idx);
            }
            TRACE_J => {
                current_state = select_j(profile, posterior_matrix, optimal_matrix, target_idx);
            }
            TRACE_E => {
                current_state = select_e(profile, optimal_matrix, target_idx, &mut profile_idx);
            }
            TRACE_B => {
                current_state = select_b(profile, optimal_matrix, target_idx);
            }
            _ => {
                panic!("bogus state in traceback")
            }
        }

        posterior_probability = get_posterior_probability(
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
            posterior_probability,
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

#[inline(always)]
pub fn get_posterior_probability(
    optimal_matrix: &DpMatrix3D,
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

#[inline(always)]
pub fn select_m(
    profile: &Profile,
    optimal_matrix: &DpMatrix3D,
    target_idx: usize,
    profile_idx: usize,
) -> usize {
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
    possible_states[argmax]
}

#[inline(always)]
pub fn select_d(
    profile: &Profile,
    optimal_matrix: &DpMatrix3D,
    target_idx: usize,
    profile_idx: usize,
) -> usize {
    let match_to_delete_path = profile
        .transition_score_delta(PROFILE_MATCH_TO_DELETE, profile_idx - 1)
        * optimal_matrix.get_match(target_idx, profile_idx - 1);
    let delete_to_delete_path = profile
        .transition_score_delta(PROFILE_DELETE_TO_DELETE, profile_idx - 1)
        * optimal_matrix.get_delete(target_idx, profile_idx - 1);

    if match_to_delete_path >= delete_to_delete_path {
        TRACE_M
    } else {
        TRACE_D
    }
}

#[inline(always)]
pub fn select_i(
    profile: &Profile,
    optimal_matrix: &DpMatrix3D,
    target_idx: usize,
    profile_idx: usize,
) -> usize {
    let match_to_insert_path = profile.transition_score_delta(PROFILE_MATCH_TO_INSERT, profile_idx)
        * optimal_matrix.get_match(target_idx - 1, profile_idx);
    let insert_to_insert_path: f32 = profile
        .transition_score_delta(PROFILE_INSERT_TO_INSERT, profile_idx)
        * optimal_matrix.get_insert(target_idx - 1, profile_idx);

    if match_to_insert_path >= insert_to_insert_path {
        TRACE_M
    } else {
        TRACE_I
    }
}

#[inline(always)]
pub fn select_n(target_idx: usize) -> usize {
    if target_idx == 0 {
        TRACE_S
    } else {
        TRACE_N
    }
}

#[inline(always)]
pub fn select_c(
    profile: &Profile,
    posterior_matrix: &DpMatrix3D,
    optimal_matrix: &DpMatrix3D,
    target_idx: usize,
) -> usize {
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

#[inline(always)]
pub fn select_j(
    profile: &Profile,
    posterior_matrix: &DpMatrix3D,
    optimal_matrix: &DpMatrix3D,
    target_idx: usize,
) -> usize {
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

#[inline(always)]
pub fn select_e(
    profile: &Profile,
    optimal_matrix: &DpMatrix3D,
    target_idx: usize,
    return_profile_idx: &mut usize,
) -> usize {
    let mut max = -f32::INFINITY;
    let mut state_max = 0;
    let mut profile_idx_max = 0;

    for profile_idx in 1..=profile.length {
        if optimal_matrix.get_match(target_idx, profile_idx) >= max {
            max = optimal_matrix.get_match(target_idx, profile_idx);
            state_max = TRACE_M;
            profile_idx_max = profile_idx;
        }
        if optimal_matrix.get_delete(target_idx, profile_idx) >= max {
            max = optimal_matrix.get_delete(target_idx, profile_idx);
            state_max = TRACE_D;
            profile_idx_max = profile_idx;
        }
    }

    *return_profile_idx = profile_idx_max;
    state_max
}

#[inline(always)]
pub fn select_b(profile: &Profile, optimal_matrix: &DpMatrix3D, target_idx: usize) -> usize {
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
