use crate::alphabet::{UTF8_DASH, UTF8_DOT, UTF8_NUMERIC, UTF8_PLUS, UTF8_SPACE};
use crate::structs::trace::constants::{TRACE_B, TRACE_D, TRACE_E, TRACE_I, TRACE_M};
use crate::structs::{Profile, Sequence, Trace};
use std::cmp::{max, min};
use std::io::Write;

use anyhow::Result;

pub struct Alignment {
    /// The name of the profile/model
    pub profile_name: String,
    /// The name of the target sequence
    pub target_name: String,
    /// The bit score of the alignment
    pub bit_score: f32,
    /// The P-value of the alignment
    pub pvalue: f64,
    /// The E-value of the alignment
    pub evalue: f64,
    /// The length of the alignment
    pub length: usize,
    /// The start coordinate of the profile (query)
    pub profile_start: usize,
    /// The end coordinate of the profile (query)
    pub profile_end: usize,
    /// The display for the profile portion of the alignment
    pub profile_string: String,
    /// The start coordinate of the target sequence
    pub target_start: usize,
    /// The end coordinate of the target sequence
    pub target_end: usize,
    /// The display for the target portion of the alignment
    pub target_string: String,
    /// The display in between the profile and target
    pub middle_string: String,
    /// The display for position-specific posterior probability bins
    pub posterior_probability_string: String,
}

enum SearchState {
    Begin,
    Alignment,
}

/// This maps a probability to a UTF8 byte (u8) to the set 0..9 or * (which represents 10)
///
/// 0.00 - 0.05 -> "0"
///
/// 0.05 - 0.15 -> "1"
///
/// ...
///
/// 0.85 - 0.95 -> "9"
///
/// 0.95 - 1.00 -> "*"
fn map_posterior_probability_to_bin_byte(probability: f32) -> u8 {
    let bin = (probability * 10.0).round();
    UTF8_NUMERIC[bin as usize]
}

/// This selects the proper character for the middle line of the alignment reporting.
fn select_middle_character(profile_byte: u8, target_byte: u8, match_emission_score: f32) -> u8 {
    if profile_byte.to_ascii_lowercase() == target_byte.to_ascii_lowercase() {
        // if we have an exact match, we just put place the matched character
        profile_byte
    } else if match_emission_score > 0.0 {
        // if we have a positive match emission score
        // (the log odds ratio), then we place a plus
        UTF8_PLUS
    } else {
        // otherwise, we just place a space
        UTF8_SPACE
    }
}

#[derive(Default)]
pub struct ScoreParams {
    pub forward_score_nats: f32,
    pub null_score_nats: f32,
    pub bias_correction_score_nats: f32,
    pub target_count: usize,
}

impl Alignment {
    pub fn from_trace(
        trace: &Trace,
        profile: &Profile,
        target: &Sequence,
        params: &ScoreParams,
    ) -> Self {
        let mut profile_bytes: Vec<u8> = vec![];
        let mut target_bytes: Vec<u8> = vec![];
        let mut mid_bytes: Vec<u8> = vec![];
        let mut posterior_probability_bytes: Vec<u8> = vec![];

        let mut profile_start: usize = 0;
        let mut target_start: usize = 0;

        let mut search_state: SearchState = SearchState::Begin;
        for trace_idx in 0..trace.length {
            let trace_state = trace.states[trace_idx];
            let profile_idx = trace.profile_idx[trace_idx];
            let profile_string_byte = profile.consensus_sequence[profile_idx];
            let target_idx = trace.target_idx[trace_idx];
            let target_string_byte = target.utf8_bytes[target_idx];

            posterior_probability_bytes.push(map_posterior_probability_to_bin_byte(
                trace.posterior_probabilities[trace_idx],
            ));

            match search_state {
                SearchState::Begin => {
                    if trace_state == TRACE_B {
                        // if we've hit a B state, then the next trace
                        // position should be the start of the alignment
                        profile_start = trace.profile_idx[trace_idx + 1];
                        target_start = trace.target_idx[trace_idx + 1];
                        search_state = SearchState::Alignment;
                    }
                }
                SearchState::Alignment => {
                    match trace_state {
                        TRACE_M => {
                            let match_emission_score = profile.match_score(
                                target.digital_bytes[target_idx] as usize,
                                profile_idx,
                            );

                            profile_bytes.push(profile_string_byte);
                            mid_bytes.push(select_middle_character(
                                profile_string_byte,
                                target_string_byte,
                                match_emission_score,
                            ));
                            target_bytes.push(target_string_byte);
                        }
                        TRACE_I => {
                            profile_bytes.push(UTF8_DOT);
                            target_bytes.push(target_string_byte);
                            mid_bytes.push(UTF8_SPACE);
                        }
                        TRACE_D => {
                            profile_bytes.push(profile.consensus_sequence[profile_idx]);
                            target_bytes.push(UTF8_DASH);
                            mid_bytes.push(UTF8_SPACE);
                        }
                        TRACE_E => {
                            // if we've hit an E state, then the previous trace
                            // position should be the end of the alignment

                            // the alignment ends at the trace position before the E node
                            let target_end = trace.target_idx[trace_idx - 1];
                            let profile_end = trace.profile_idx[trace_idx - 1];

                            let aligned_target_length = target_end - target_start + 1;
                            let unaligned_target_length =
                                (target.length - aligned_target_length) as f32;

                            // sum up the loop transitions to the N and/or C states
                            // once for every position in the target sequence that
                            // isn't accounted for by the alignment that we produced
                            let n_and_c_state_correction_nats = unaligned_target_length
                                * (target.length as f32 / (target.length as f32 + 3.0)).ln();

                            // everything has been in nats up to here
                            let nat_score = params.forward_score_nats
                                + n_and_c_state_correction_nats
                                - (params.null_score_nats + params.bias_correction_score_nats);

                            // convert to bits by dividing by ln(2)
                            let bit_score = nat_score / std::f32::consts::LN_2;

                            // TODO: double check these calculations
                            let pvalue = (-profile.forward_lambda as f64
                                * (bit_score as f64 - profile.forward_tau as f64))
                                .exp();

                            let evalue = pvalue * params.target_count as f64;

                            return Alignment {
                                profile_name: profile.name.clone(),
                                target_name: target.name.clone(),
                                bit_score,
                                pvalue,
                                evalue,
                                length: profile_bytes.len(),
                                profile_start,
                                profile_end,
                                profile_string: String::from_utf8(profile_bytes).unwrap(),
                                target_start,
                                target_end,
                                target_string: String::from_utf8(target_bytes).unwrap(),
                                middle_string: String::from_utf8(mid_bytes).unwrap(),
                                posterior_probability_string: String::from_utf8(
                                    posterior_probability_bytes,
                                )
                                .unwrap(),
                            };
                        }
                        _ => {
                            // TODO: Error
                            panic!("unknown state in Alignment::from_trace()")
                        }
                    }
                }
            }
        }
        // TODO: Error
        panic!("failed to produce an Alignment in Alignment::from_trace()");
    }

    pub fn dump(&self, out: &mut impl Write) -> Result<()> {
        let mut start_offset: usize = 0;
        let mut end_offset: usize = 80;

        let name_width = max(self.profile_name.len(), self.target_name.len());

        // write the score line
        writeln!(
            out,
            "==  score: {:3.1} bits;  E-value: {:1.1e}",
            self.bit_score, self.evalue
        )?;

        while start_offset <= self.length {
            start_offset = min(start_offset, self.length);
            end_offset = min(end_offset, self.length);

            // TODO: consensus structure line?

            // write the profile sequence
            writeln!(
                out,
                "{:>W$} {:5} {} {:<5}",
                self.profile_name,
                self.profile_start + start_offset,
                &self.profile_string[start_offset..end_offset],
                self.profile_start + end_offset - 1,
                W = name_width
            )?;

            // write the middle line
            writeln!(
                out,
                "{:W$} {:5} {}",
                "",
                "",
                &self.middle_string[start_offset..end_offset],
                W = name_width
            )?;

            // write the target sequence
            writeln!(
                out,
                "{:>W$} {:5} {} {:<5}",
                self.target_name,
                self.target_start + start_offset,
                &self.target_string[start_offset..end_offset],
                self.target_start + end_offset - 1,
                W = name_width
            )?;

            // write the position-specific posterior probabilities
            writeln!(
                out,
                "{:W$} {:5} {}",
                "",
                "",
                &self.posterior_probability_string[start_offset..end_offset],
                W = name_width
            )?;

            writeln!(out)?;

            start_offset += 80;
            end_offset += 80;
        }

        Ok(())
    }
}
