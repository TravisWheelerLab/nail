use crate::align::{e_value, p_value, Bits, Score};
use crate::alphabet::{UTF8_DASH, UTF8_DOT, UTF8_NUMERIC, UTF8_PLUS, UTF8_SPACE};
use crate::output::output_tabular::{Field, TableFormat};
use crate::structs::{Profile, Sequence};
use std::cmp::{max, min};

use super::trace::TraceStep;
use super::Trace;

pub struct Alignment {
    /// The length of the alignment
    pub length: usize,
    /// The start coordinate of the profile (query)
    pub profile_start: usize,
    /// The end coordinate of the profile (query)
    pub profile_end: usize,
    /// The start coordinate of the target sequence
    pub target_start: usize,
    /// The end coordinate of the target sequence
    pub target_end: usize,

    // optional fields
    /// The name of the profile/model
    pub profile_name: Option<String>,
    /// The name of the target sequence
    pub target_name: Option<String>,
    /// The final bit score of the alignment (after all bias adjustments)
    pub score_bits: Option<Bits>,
    /// The raw bit score of the alignment (before any bias adjustments)
    pub raw_score_bits: Option<Bits>,
    /// The bit score of the length bias adjustment
    pub null_one: Option<Bits>,
    /// The bit score of the composition bias adjustment
    pub null_two: Option<Bits>,
    /// The P-value of the alignment
    pub p_value: Option<f64>,
    /// The E-value of the alignment
    pub e_value: Option<f64>,
    /// The fraction of dynamic programming cells filled in during alignment.
    pub cell_fraction: Option<f32>,

    // display strings
    /// The display for the profile portion of the alignment
    pub profile_string: Option<String>,
    /// The display for the target portion of the alignment
    pub target_string: Option<String>,
    /// The display in between the profile and target
    pub middle_string: Option<String>,
    /// The display for position-specific posterior probability bins
    pub posterior_string: Option<String>,
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

#[derive(Default)]
pub struct AlignmentBuilder<'a> {
    trace: Vec<TraceStep>,
    profile: Option<&'a Profile>,
    target: Option<&'a Sequence>,
    target_count: Option<usize>,
    forward_score: Option<Bits>,
    null_one: Option<Bits>,
    null_two: Option<Bits>,
    cell_fraction: Option<f32>,
}

impl<'a> AlignmentBuilder<'a> {
    pub fn new(trace: &Trace) -> Self {
        let trace = trace
            .iter()
            .filter(|s| {
                s.state == Trace::M_STATE || s.state == Trace::I_STATE || s.state == Trace::D_STATE
            })
            .collect();

        Self {
            trace,
            ..Default::default()
        }
    }

    pub fn with_profile(mut self, profile: &'a Profile) -> Self {
        self.profile = Some(profile);
        self
    }

    pub fn with_target(mut self, target: &'a Sequence) -> Self {
        self.target = Some(target);
        self
    }

    pub fn with_target_count(mut self, count: usize) -> Self {
        self.target_count = Some(count);
        self
    }

    pub fn with_forward_score(mut self, score: impl Score) -> Self {
        self.forward_score = Some(score.bits());
        self
    }

    pub fn with_null_one(mut self, score: impl Score) -> Self {
        self.null_one = Some(score.bits());
        self
    }

    pub fn with_null_two(mut self, score: impl Score) -> Self {
        self.null_two = Some(score.bits());
        self
    }

    pub fn with_cell_fraction(mut self, cell_fraction: f32) -> Self {
        self.cell_fraction = Some(cell_fraction);
        self
    }

    pub fn build(self) -> anyhow::Result<Alignment> {
        let length = self.trace.len();

        let first = self.trace.first().unwrap();
        let last = self.trace.last().unwrap();

        let profile_start = first.profile_idx;
        let profile_end = last.profile_idx;
        let target_start = first.target_idx;
        let target_end = last.target_idx;

        let score = match self.forward_score {
            Some(mut score) => {
                // TODO: decide if we actually want to do this correction
                if let Some(target) = self.target {
                    let aligned_target_length = target_end - target_start + 1;
                    let unaligned_target_length = (target.length - aligned_target_length) as f32;

                    // sum up the loop transitions to the N and/or C states
                    // once for every position in the target sequence that
                    // isn't accounted for by the alignment that we produced
                    let correction = Bits(
                        unaligned_target_length
                            * (target.length as f32 / (target.length as f32 + 3.0)).log2(),
                    );

                    score = score + correction
                }

                if let Some(null_one) = self.null_one {
                    score = score - null_one
                }

                if let Some(null_two) = self.null_two {
                    score = score - null_two
                }

                Some(score)
            }
            None => None,
        };

        let p_value = match (score, self.profile) {
            (Some(score), Some(profile)) => {
                Some(p_value(score, profile.forward_lambda, profile.forward_tau))
            }
            (_, _) => None,
        };

        let e_value = p_value.map(|p_value| e_value(p_value, self.target_count.unwrap_or(1)));

        let profile_name = self.profile.map(|profile| profile.name.clone());

        let target_name = self.target.map(|target| target.name.clone());

        let (profile_string, target_string, middle_string, posterior_string) =
            match (self.profile, self.target) {
                (Some(profile), Some(target)) => {
                    let mut profile_bytes = vec![];
                    let mut target_bytes = vec![];
                    let mut middle_bytes = vec![];
                    let mut posteriors = vec![];

                    self.trace.iter().for_each(|step| {
                        posteriors.push(map_posterior_probability_to_bin_byte(step.posterior));
                        let profile_byte = profile.consensus_sequence_bytes_utf8[step.profile_idx];
                        let target_byte = target.utf8_bytes[step.target_idx];

                        match step.state {
                            Trace::I_STATE => {
                                profile_bytes.push(UTF8_DOT);
                                target_bytes.push(target_byte);
                                middle_bytes.push(UTF8_SPACE);
                            }
                            Trace::D_STATE => {
                                profile_bytes.push(profile_byte);
                                target_bytes.push(UTF8_DASH);
                                middle_bytes.push(UTF8_SPACE);
                            }
                            Trace::M_STATE => {
                                let target_byte_digital = target.digital_bytes[step.target_idx];

                                profile_bytes.push(profile_byte);
                                target_bytes.push(target_byte);

                                if profile_byte == target_byte {
                                    middle_bytes.push(profile_byte);
                                } else if profile
                                    .match_score(target_byte_digital as usize, step.profile_idx)
                                    > 0.0
                                {
                                    middle_bytes.push(UTF8_PLUS);
                                } else {
                                    middle_bytes.push(UTF8_SPACE);
                                }
                            }
                            _ => {
                                panic!("invalid trace state in AlignmentBuilder: {}", step.state)
                            }
                        }
                    });

                    let profile_string = String::from_utf8(profile_bytes)?;
                    let target_string = String::from_utf8(target_bytes)?;
                    let middle_string = String::from_utf8(middle_bytes)?;
                    let posterior_string = String::from_utf8(posteriors)?;
                    (
                        Some(profile_string),
                        Some(target_string),
                        Some(middle_string),
                        Some(posterior_string),
                    )
                }
                _ => (None, None, None, None),
            };

        Ok(Alignment {
            profile_name,
            target_name,
            score_bits: score,
            raw_score_bits: self.forward_score,
            null_one: self.null_one,
            null_two: self.null_two,
            p_value,
            e_value,
            length,
            cell_fraction: self.cell_fraction,
            profile_start,
            profile_end,
            target_start,
            target_end,
            profile_string,
            target_string,
            middle_string,
            posterior_string,
        })
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
            let profile_idx = trace.profile_indices[trace_idx];
            let profile_string_byte = profile.consensus_sequence_bytes_utf8[profile_idx];
            let target_idx = trace.target_indices[trace_idx];
            let target_string_byte = target.utf8_bytes[target_idx];

            posterior_probability_bytes.push(map_posterior_probability_to_bin_byte(
                trace.posterior_probabilities[trace_idx],
            ));

            match search_state {
                SearchState::Begin => {
                    if trace_state == Trace::B_STATE {
                        // if we've hit a B state, then the next trace
                        // position should be the start of the alignment
                        profile_start = trace.profile_indices[trace_idx + 1];
                        target_start = trace.target_indices[trace_idx + 1];
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
                            profile_bytes.push(profile.consensus_sequence_bytes_utf8[profile_idx]);
                            target_bytes.push(UTF8_DASH);
                            mid_bytes.push(UTF8_SPACE);
                        }
                        Trace::E_STATE => {
                            // if we've hit an E state, then the previous trace
                            // position should be the end of the alignment

                            // the alignment ends at the trace position before the E node
                            let target_end = trace.target_indices[trace_idx - 1];
                            let profile_end = trace.profile_indices[trace_idx - 1];

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
                                profile_name: Some(profile.name.clone()),
                                target_name: Some(target.name.clone()),
                                score_bits: Some(Bits(score_bits)),
                                raw_score_bits: Some(Bits(raw_score_bits)),
                                null_one: Some(Bits(length_bias_bits)),
                                null_two: Some(Bits(composition_bias_bits)),
                                p_value: Some(pvalue),
                                e_value: Some(evalue),
                                length: profile_bytes.len(),
                                profile_start,
                                profile_end,
                                target_start,
                                target_end,
                                cell_fraction: None,
                                profile_string: Some(String::from_utf8(profile_bytes).unwrap()),
                                target_string: Some(String::from_utf8(target_bytes).unwrap()),
                                middle_string: Some(String::from_utf8(mid_bytes).unwrap()),
                                posterior_string: Some(
                                    String::from_utf8(posterior_probability_bytes).unwrap(),
                                ),
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

    pub fn tab_string_formatted(&self, format: &TableFormat) -> String {
        let mut tab_string = String::new();

        format
            .fields
            .iter()
            .zip(format.widths.iter())
            .for_each(|(field, width)| {
                let val = field.extract_from(self);
                tab_string = format!("{tab_string}{val:width$} ", width = width)
            });

        // remove the last space
        tab_string.pop();

        tab_string
    }

    pub fn ali_string(&self) -> String {
        let (profile_string, target_string, middle_string, posterior_string) = match (
            &self.profile_string,
            &self.target_string,
            &self.middle_string,
            &self.posterior_string,
        ) {
            (Some(a), Some(b), Some(c), Some(d)) => (a, b, c, d),
            _ => panic!(),
        };

        let mut ali_string = String::new();
        let mut start_offset: usize = 0;
        let mut end_offset: usize = 80;

        let name_width = max(
            Field::Query.extract_from(self).len(),
            Field::Target.extract_from(self).len(),
        );

        // score line
        ali_string.push_str(&format!(
            "==  score: {} bits;  E-value: {}\n",
            Field::Score.extract_from(self),
            Field::Evalue.extract_from(self)
        ));

        while start_offset <= self.length {
            start_offset = min(start_offset, self.length);
            end_offset = min(end_offset, self.length);

            // profile sequence
            ali_string.push_str(&format!(
                "{:>W$} {:5} {} {:<5}\n",
                Field::Query.extract_from(self).len(),
                self.profile_start + start_offset,
                &profile_string[start_offset..end_offset],
                self.profile_start + end_offset - 1,
                W = name_width
            ));

            // middle line
            ali_string.push_str(&format!(
                "{:W$} {:5} {}\n",
                "",
                "",
                &middle_string[start_offset..end_offset],
                W = name_width
            ));

            // target sequence
            ali_string.push_str(&format!(
                "{:>W$} {:5} {} {:<5}\n",
                Field::Target.extract_from(self).len(),
                self.target_start + start_offset,
                &target_string[start_offset..end_offset],
                self.target_start + end_offset - 1,
                W = name_width
            ));

            // position-specific posterior probabilities
            ali_string.push_str(&format!(
                "{:W$} {:5} {}\n\n",
                "",
                "",
                &posterior_string[start_offset..end_offset],
                W = name_width
            ));

            start_offset += 80;
            end_offset += 80;
        }

        ali_string
    }
}
