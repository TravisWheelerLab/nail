use crate::alphabet::{UTF8_DASH, UTF8_DOT, UTF8_NUMERIC, UTF8_PLUS, UTF8_SPACE};
use crate::output::output_tabular::{Field, TableFormat};
use crate::structs::{Profile, Sequence};
use std::cmp::{max, min};

use super::Trace;

pub struct Alignment {
    /// The name of the profile/model
    pub profile_name: String,
    /// The name of the target sequence
    pub target_name: String,
    /// The final bit score of the alignment (after all bias adjustments)
    pub score_bits: f32,
    /// The raw bit score of the alignment (before any bias adjustments)
    pub raw_score_bits: f32,
    /// The bit score of the length bias adjustment
    pub length_bias_bits: f32,
    /// The bit score of the composition bias adjustment
    pub composition_bias_bits: f32,
    /// The P-value of the alignment
    pub pvalue: f64,
    /// The E-value of the alignment
    pub evalue: f64,
    /// The length of the alignment
    pub length: usize,
    /// The fraction of dynamic programming cells filled in during alignment.
    pub cell_fraction: Option<f32>,
    /// The start coordinate of the profile (query)
    pub profile_start: usize,
    /// The end coordinate of the profile (query)
    pub profile_end: usize,
    /// The start coordinate of the target sequence
    pub target_start: usize,
    /// The end coordinate of the target sequence
    pub target_end: usize,

    // display strings
    // ---------------
    /// The display for the profile portion of the alignment
    pub profile_string: String,
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

#[derive(Default, Clone)]
pub struct ScoreParams {
    pub forward_score_nats: f32,
    pub length_bias_score_nats: f32,
    pub composition_bias_score_nats: f32,
    pub target_count: usize,
}

impl ScoreParams {
    pub fn new(target_count: usize) -> Self {
        Self {
            forward_score_nats: 0.0,
            length_bias_score_nats: 0.0,
            composition_bias_score_nats: 0.0,
            target_count,
        }
    }
}

struct AlignmentBuilder {
    profile_indices: Vec<usize>,
    target_indices: Vec<usize>,
    posteriors: Vec<f32>,
}

impl AlignmentBuilder {
    fn new(trace: &Trace) -> Self {
        let mut profile_indices = vec![];
        let mut target_indices = vec![];
        let mut posteriors = vec![];

        trace
            .iter()
            .filter(|s| {
                s.state == Trace::M_STATE || s.state == Trace::I_STATE || s.state == Trace::D_STATE
            })
            .for_each(|s| {
                profile_indices.push(s.profile_idx);
                target_indices.push(s.target_idx);
                posteriors.push(s.posterior_probability);
            });

        Self {
            profile_indices,
            target_indices,
            posteriors,
        }
    }

    fn build(self) -> Alignment {
        Alignment {
            profile_name: todo!(),
            target_name: todo!(),
            score_bits: todo!(),
            raw_score_bits: todo!(),
            length_bias_bits: todo!(),
            composition_bias_bits: todo!(),
            pvalue: todo!(),
            evalue: todo!(),
            length: todo!(),
            cell_fraction: todo!(),
            profile_start: todo!(),
            profile_end: todo!(),
            target_start: todo!(),
            target_end: todo!(),
            profile_string: todo!(),
            target_string: todo!(),
            middle_string: todo!(),
            posterior_probability_string: todo!(),
        }
    }
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
                    if trace_state == Trace::B_STATE {
                        // if we've hit a B state, then the next trace
                        // position should be the start of the alignment
                        profile_start = trace.profile_idx[trace_idx + 1];
                        target_start = trace.target_idx[trace_idx + 1];
                        search_state = SearchState::Alignment;
                    }
                }
                SearchState::Alignment => {
                    match trace_state {
                        Trace::M_STATE => {
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
                        Trace::I_STATE => {
                            profile_bytes.push(UTF8_DOT);
                            target_bytes.push(target_string_byte);
                            mid_bytes.push(UTF8_SPACE);
                        }
                        Trace::D_STATE => {
                            profile_bytes.push(profile.consensus_sequence[profile_idx]);
                            target_bytes.push(UTF8_DASH);
                            mid_bytes.push(UTF8_SPACE);
                        }
                        Trace::E_STATE => {
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

                            let raw_score_bits = (params.forward_score_nats
                                + n_and_c_state_correction_nats)
                                / std::f32::consts::LN_2;

                            let length_bias_bits =
                                params.length_bias_score_nats / std::f32::consts::LN_2;

                            let composition_bias_bits =
                                params.composition_bias_score_nats / std::f32::consts::LN_2;

                            let score_bits =
                                raw_score_bits - length_bias_bits - composition_bias_bits;

                            // TODO: double check these calculations
                            let pvalue = (-profile.forward_lambda as f64
                                * (score_bits as f64 - profile.forward_tau as f64))
                                .exp();

                            let evalue = pvalue * params.target_count as f64;

                            return Alignment {
                                profile_name: profile.name.clone(),
                                target_name: target.name.clone(),
                                score_bits,
                                raw_score_bits,
                                length_bias_bits,
                                composition_bias_bits,
                                pvalue,
                                evalue,
                                length: profile_bytes.len(),
                                profile_start,
                                profile_end,
                                target_start,
                                target_end,
                                cell_fraction: None,
                                profile_string: String::from_utf8(profile_bytes).unwrap(),
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

    pub const TAB_HEADER: &'static str = "#target\tquery\ttarget start\ttarget end\tprofile start\tprofile end\tscore\tcomposition bias\tE-value\tcell fraction";

    pub fn tab_string(&self) -> String {
        format!(
            "{} {} {} {} {} {} {:.2} {:.2} {:.1e} {}",
            self.target_name,
            self.profile_name,
            self.target_start,
            self.target_end,
            self.profile_start,
            self.profile_end,
            self.score_bits,
            self.composition_bias_bits,
            self.evalue,
            match self.cell_fraction {
                Some(frac) => format!("{:.1e}", frac),
                None => "-".to_string(),
            },
        )
    }

    pub fn tab_string_formatted(&self, format: &TableFormat) -> String {
        let mut tab_string = String::new();

        format
            .fields
            .iter()
            .zip(format.widths.iter())
            .for_each(|(field, width)| {
                let val = match field {
                    Field::Target => self.target_name.clone(),
                    Field::Query => self.profile_name.clone(),
                    Field::TargetStart => self.target_start.to_string(),
                    Field::TargetEnd => self.target_end.to_string(),
                    Field::QueryStart => self.profile_start.to_string(),
                    Field::QueryEnd => self.profile_end.to_string(),
                    Field::Score => format!("{:.2}", self.score_bits),
                    Field::CompBias => format!("{:.2}", self.composition_bias_bits),
                    Field::Evalue => format!("{:.1e}", self.evalue),
                    Field::CellFrac => match self.cell_fraction {
                        Some(frac) => format!("{:.1e}", frac),
                        None => "-".to_string(),
                    },
                };
                tab_string = format!("{tab_string}{val:width$} ", width = width)
            });

        // remove the last space
        tab_string.pop();

        tab_string
    }

    pub fn ali_string(&self) -> String {
        let mut ali_string = String::new();
        let mut start_offset: usize = 0;
        let mut end_offset: usize = 80;

        let name_width = max(self.profile_name.len(), self.target_name.len());

        // score line
        ali_string.push_str(&format!(
            "==  score: {:3.1} bits;  E-value: {:1.1e}\n",
            self.score_bits, self.evalue
        ));

        while start_offset <= self.length {
            start_offset = min(start_offset, self.length);
            end_offset = min(end_offset, self.length);

            // profile sequence
            ali_string.push_str(&format!(
                "{:>W$} {:5} {} {:<5}\n",
                self.profile_name,
                self.profile_start + start_offset,
                &self.profile_string[start_offset..end_offset],
                self.profile_start + end_offset - 1,
                W = name_width
            ));

            // middle line
            ali_string.push_str(&format!(
                "{:W$} {:5} {}\n",
                "",
                "",
                &self.middle_string[start_offset..end_offset],
                W = name_width
            ));

            // target sequence
            ali_string.push_str(&format!(
                "{:>W$} {:5} {} {:<5}\n",
                self.target_name,
                self.target_start + start_offset,
                &self.target_string[start_offset..end_offset],
                self.target_start + end_offset - 1,
                W = name_width
            ));

            // position-specific posterior probabilities
            ali_string.push_str(&format!(
                "{:W$} {:5} {}\n\n",
                "",
                "",
                &self.posterior_probability_string[start_offset..end_offset],
                W = name_width
            ));

            start_offset += 80;
            end_offset += 80;
        }

        ali_string
    }
}
