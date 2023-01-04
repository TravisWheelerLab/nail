use crate::alphabet::{string_from_amino_bytes, DASH_UTF8, DOT_UTF8, PIPE_UTF8, SPACE_UTF8};
use crate::structs::trace::constants::{
    TRACE_B, TRACE_C, TRACE_D, TRACE_E, TRACE_I, TRACE_M, TRACE_N, TRACE_S, TRACE_T,
};
use crate::structs::{Profile, Sequence, Trace};

pub struct Alignment {
    // pub profile_name: String,
    // pub target_name: String,
    // pub score: f32,
    pub profile_start: usize,
    pub profile_end: usize,
    pub profile_string: String,
    pub target_start: usize,
    pub target_end: usize,
    pub target_string: String,
    pub middle_string: String,
}

enum SearchState {
    TraceBegin,
    AlignmentBegin,
    AlignmentCore,
    TraceEnd,
}

impl Alignment {
    pub fn new(trace: &Trace, profile: &Profile, target: &Sequence) -> Self {
        let mut profile_bytes: Vec<u8> = vec![];
        let mut target_bytes: Vec<u8> = vec![];
        let mut mid_bytes: Vec<u8> = vec![];

        let mut profile_start: usize = 0;
        let mut profile_end: usize = 0;
        let mut target_start: usize = 0;
        let mut target_end: usize = 0;

        let mut search_state: SearchState = SearchState::TraceBegin;

        for trace_idx in 0..trace.length {
            let trace_state = trace.states[trace_idx];
            let profile_idx = trace.profile_idx[trace_idx];
            let target_idx = trace.target_idx[trace_idx];
            match search_state {
                SearchState::TraceBegin => match trace_state {
                    TRACE_S | TRACE_N => {}
                    TRACE_B => {
                        search_state = SearchState::AlignmentBegin;
                    }
                    _ => {
                        panic!();
                    }
                },
                SearchState::AlignmentBegin => match trace_state {
                    TRACE_M => {
                        profile_start = profile_idx;
                        target_start = target_idx;
                        profile_bytes.push(profile.consensus_sequence[profile_idx]);
                        target_bytes.push(target.data[target_idx]);
                        mid_bytes.push(PIPE_UTF8);
                        search_state = SearchState::AlignmentCore;
                    }
                    TRACE_I => {
                        profile_start = profile_idx;
                        target_start = target_idx;
                        profile_bytes.push(DOT_UTF8);
                        mid_bytes.push(SPACE_UTF8);
                        target_bytes.push(target.data[target_idx]);
                        search_state = SearchState::AlignmentCore;
                    }
                    TRACE_D => {
                        profile_start = profile_idx;
                        target_start = target_idx;
                        profile_bytes.push(profile.consensus_sequence[profile_idx]);
                        target_bytes.push(DASH_UTF8);
                        search_state = SearchState::AlignmentCore;
                    }
                    _ => {
                        panic!()
                    }
                },
                SearchState::AlignmentCore => match trace_state {
                    TRACE_M => {
                        profile_bytes.push(profile.consensus_sequence[profile_idx]);
                        target_bytes.push(target.data[target_idx]);
                        mid_bytes.push(PIPE_UTF8);
                    }
                    TRACE_I => {
                        profile_bytes.push(DOT_UTF8);
                        target_bytes.push(target.data[target_idx]);
                        mid_bytes.push(SPACE_UTF8);
                    }
                    TRACE_D => {
                        profile_bytes.push(profile.consensus_sequence[profile_idx]);
                        target_bytes.push(DASH_UTF8);
                        mid_bytes.push(SPACE_UTF8);
                    }
                    TRACE_E => {
                        profile_end = trace.profile_idx[trace_idx - 1];
                        target_end = trace.target_idx[trace_idx - 1];
                        search_state = SearchState::TraceEnd;
                    }
                    _ => {
                        panic!()
                    }
                },
                SearchState::TraceEnd => {}
            }
        }

        Alignment {
            profile_start,
            profile_end,
            profile_string: String::from_utf8(profile_bytes).unwrap(),
            target_start,
            target_end,
            target_string: string_from_amino_bytes(&target_bytes),
            middle_string: String::from_utf8(mid_bytes).unwrap(),
        }
    }
}
