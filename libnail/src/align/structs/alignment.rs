use anyhow::bail;

use crate::align::{e_value, p_value, Bits, Score};
use crate::alphabet::{UTF8_DASH, UTF8_DOT, UTF8_NUMERIC, UTF8_PLUS, UTF8_SPACE};
use crate::output::output_tabular::{Field, TableFormat};
use crate::structs::{Profile, Sequence};
use std::cmp::{max, min};

use super::Trace;

#[derive(Default)]
pub struct Boundaries {
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
}

pub struct Scores {
    /// The Forward score (no bias adjustment)
    pub forward_score: Bits,
    /// The Forward P-value (no bias adjustment)
    pub forward_p_value: f64,
    /// The final bit score of the alignment (after bias adjustments)
    pub bit_score: Bits,
    /// The bit score of the composition bias adjustment
    pub null_two_score: Option<Bits>,
    /// The P-value of the alignment
    pub p_value: f64,
    /// The E-value of the alignment
    pub e_value: f64,
}

pub struct CellStats {
    /// The number of dynamic programming cells filled in during alignment.
    pub count: usize,
    /// The fraction of dynamic programming cells filled in during alignment.
    pub fraction: f32,
}

pub struct DisplayStrings {
    /// The display for the profile portion of the alignment
    pub profile_string: String,
    /// The display for the target portion of the alignment
    pub target_string: String,
    /// The display in between the profile and target
    pub middle_string: String,
    /// The display for position-specific posterior probability bins
    pub posterior_string: String,
}

pub struct Alignment {
    /// The name of the profile/model
    pub profile_name: Option<String>,
    /// The name of the target sequence
    pub target_name: Option<String>,
    /// The boundaries of the alignment
    pub boundaries: Option<Boundaries>,
    /// The bitscores and significance metrics of the alignment
    pub scores: Scores,
    /// The metrics for the sparse dynamic programming matrix used to compute the alignment
    pub cell_stats: Option<CellStats>,
    /// The strings used for alignment display
    pub display_strings: Option<DisplayStrings>,
}

impl AsRef<Alignment> for &Alignment {
    fn as_ref(&self) -> &Alignment {
        self
    }
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
    target: Option<&'a Sequence>,
    profile: Option<&'a Profile>,
    trace: Option<&'a Trace>,
    database_size: Option<usize>,
    forward_score: Option<Bits>,
    null_two: Option<Bits>,
    cell_count: Option<usize>,
}

impl<'a> AlignmentBuilder<'a> {
    pub fn with_trace(mut self, trace: &'a Trace) -> Self {
        self.trace = Some(trace);
        self
    }

    pub fn with_profile(mut self, profile: &'a Profile) -> Self {
        self.profile = Some(profile);
        self
    }

    pub fn with_target(mut self, target: &'a Sequence) -> Self {
        self.target = Some(target);
        self
    }

    pub fn with_database_size(mut self, count: usize) -> Self {
        self.database_size = Some(count);
        self
    }

    pub fn with_forward_score(mut self, score: impl Score) -> Self {
        self.forward_score = Some(score.bits());
        self
    }

    pub fn with_null_two(mut self, score: Option<impl Score>) -> Self {
        // we don't allow positive bias composition
        self.null_two = score.map(|s| s.bits().min(Bits(0.0)));
        self
    }

    pub fn with_cell_count(mut self, cell_count: usize) -> Self {
        self.cell_count = Some(cell_count);
        self
    }

    pub fn build(self) -> anyhow::Result<Alignment> {
        let scores = match self.forward_score {
            Some(forward_score) => {
                let bit_score = match self.null_two {
                    Some(null_two) => forward_score - null_two,
                    None => forward_score,
                };

                let (forward_p_value, p_value) = match self.profile {
                    Some(profile) => (
                        p_value(forward_score, profile.forward_lambda, profile.forward_tau),
                        p_value(bit_score, profile.forward_lambda, profile.forward_tau),
                    ),
                    _ => bail!("Profile missing during Alignment construction"),
                };

                // default to a database size of 1
                let e_value = e_value(p_value, self.database_size.unwrap_or(1));

                Scores {
                    forward_score,
                    forward_p_value,
                    bit_score,
                    null_two_score: self.null_two,
                    p_value,
                    e_value,
                }
            }
            None => bail!("score missing during Alignment construction"),
        };

        let boundaries = match self.trace {
            Some(trace) => {
                let length = trace.core_len();
                let first = trace.first_core();
                let last = trace.last_core();

                match (first, last) {
                    (Some(first), Some(last)) => Some(Boundaries {
                        length,
                        profile_start: first.profile_idx,
                        profile_end: last.profile_idx,
                        target_start: first.target_idx,
                        target_end: last.target_idx,
                    }),
                    // if we have a trace, but it has no core
                    // steps, set all Boundary fields to 0
                    (_, _) => Some(Boundaries::default()),
                }
            }
            None => None,
        };

        let cell_stats = match self.cell_count {
            Some(count) => {
                let fraction = match (self.target, self.profile) {
                    (Some(target), Some(profile)) => {
                        count as f32 / (target.length * profile.length) as f32
                    }
                    (_, _) => 1.0,
                };
                Some(CellStats { count, fraction })
            }
            None => None,
        };

        let display_strings = match (self.trace, self.profile, self.target) {
            (Some(trace), Some(profile), Some(target)) => {
                let mut profile_bytes = vec![];
                let mut target_bytes = vec![];
                let mut middle_bytes = vec![];
                let mut posteriors = vec![];

                trace
                    .iter()
                    .filter(|s| {
                        s.state == Trace::M_STATE
                            || s.state == Trace::I_STATE
                            || s.state == Trace::D_STATE
                    })
                    .for_each(|step| {
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

                                if profile_byte.to_ascii_lowercase()
                                    == target_byte.to_ascii_lowercase()
                                {
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

                Some(DisplayStrings {
                    profile_string,
                    target_string,
                    middle_string,
                    posterior_string,
                })
            }
            _ => None,
        };

        Ok(Alignment {
            profile_name: self.profile.map(|profile| profile.name.clone()),
            target_name: self.target.map(|target| target.name.clone()),
            boundaries,
            scores,
            cell_stats,
            display_strings,
        })
    }
}

impl Alignment {
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
        match (&self.display_strings, &self.boundaries) {
            (Some(display), Some(boundaries)) => {
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

                while start_offset <= boundaries.length {
                    start_offset = min(start_offset, boundaries.length);
                    end_offset = min(end_offset, boundaries.length);

                    // profile sequence
                    ali_string.push_str(&format!(
                        "{:>W$} {:5} {} {:<5}\n",
                        Field::Query.extract_from(self),
                        boundaries.profile_start + start_offset,
                        &display.profile_string[start_offset..end_offset],
                        boundaries.profile_start + end_offset - 1,
                        W = name_width
                    ));

                    // middle line
                    ali_string.push_str(&format!(
                        "{:W$} {:5} {}\n",
                        "",
                        "",
                        &display.middle_string[start_offset..end_offset],
                        W = name_width
                    ));

                    // target sequence
                    ali_string.push_str(&format!(
                        "{:>W$} {:5} {} {:<5}\n",
                        Field::Target.extract_from(self),
                        boundaries.target_start + start_offset,
                        &display.target_string[start_offset..end_offset],
                        boundaries.target_start + end_offset - 1,
                        W = name_width
                    ));

                    // position-specific posterior probabilities
                    ali_string.push_str(&format!(
                        "{:W$} {:5} {}\n\n",
                        "",
                        "",
                        &display.posterior_string[start_offset..end_offset],
                        W = name_width
                    ));

                    start_offset += 80;
                    end_offset += 80;
                }

                ali_string
            }
            (_, _) => panic!(),
        }
    }
}
