use std::fmt::Display;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use crate::structs::profile::constants::MAX_DEGENERATE_ALPHABET_SIZE;
use anyhow::{Context, Result};
use lazy_static::lazy_static;
use regex::Regex;
use thiserror::Error;

pub mod constants {
    // these constants describe indices of transitions
    pub const HMM_MATCH_TO_MATCH: usize = 0;
    pub const HMM_MATCH_TO_INSERT: usize = 1;
    pub const HMM_MATCH_TO_DELETE: usize = 2;
    pub const HMM_INSERT_TO_MATCH: usize = 3;
    pub const HMM_INSERT_TO_INSERT: usize = 4;
    pub const HMM_DELETE_TO_MATCH: usize = 5;
    pub const HMM_DELETE_TO_DELETE: usize = 6;
}

// local constants for parsing flags
const P7_HEADER_FORMAT_FLAG: &'static str = "HMMER3/f";
const P7_HEADER_NAME_FLAG: &'static str = "NAME";
const P7_HEADER_ACCESSION_FLAG: &'static str = "ACC";
const P7_HEADER_DESCRIPTION_FLAG: &'static str = "DESC";
const P7_HEADER_LENGTH_FLAG: &'static str = "LENG";
const P7_HEADER_MAXL_FLAG: &'static str = "MAXL";
const P7_HEADER_ALPHABET_FLAG: &'static str = "ALPH";
const P7_HEADER_REFERENCE_FLAG: &'static str = "RF";
const P7_HEADER_MASK_FLAG: &'static str = "MM";
const P7_HEADER_CONSENSUS_RESIDUE_FLAG: &'static str = "CONS";
const P7_HEADER_CONSENSUS_STRUCTURE_FLAG: &'static str = "CS";
const P7_HEADER_MAP_FLAG: &'static str = "MAP";
const P7_HEADER_DATE_FLAG: &'static str = "DATE";
const P7_HEADER_COMMAND_FLAG: &'static str = "COM";
const P7_HEADER_NSEQ_FLAG: &'static str = "NSEQ";
const P7_HEADER_EFFN_FLAG: &'static str = "EFFN";
const P7_HEADER_CHECKSUM_FLAG: &'static str = "CKSUM";
const P7_HEADER_GATHERING_FLAG: &'static str = "GA";
const P7_HEADER_TRUSTED_FLAG: &'static str = "TC";
const P7_HEADER_NOISE_FLAG: &'static str = "NC";
const P7_HEADER_STATS_FLAG: &'static str = "STATS";
const P7_HEADER_STATS_MSV_FLAG: &'static str = "MSV";
const P7_HEADER_STATS_VITERBI_FLAG: &'static str = "VITERBI";
const P7_HEADER_STATS_FORWARD_FLAG: &'static str = "FORWARD";

const P7_BODY_HMM_MODEL_START_FLAG: &'static str = "HMM";
const P7_BODY_COMPO_FLAG: &'static str = "COMPO";
const P7_BODY_END_FLAG: &'static str = "//";

// this static regex is used to find float strings
lazy_static! {
    static ref FLOAT_RE: Regex = Regex::new(r"\d\.*\d*").unwrap();
}

enum ParserState {
    Idle,
    Header,
    ModelHead,
    ModelBody,
}

enum ModelParserState {
    MatchEmissions,
    InsertEmissions,
    StateTransitions,
}

/// An Error that is thrown when an unknown header is found
/// when parsing the header of a P7HMM
#[derive(Error, Debug)]
#[error("unknown header flag")]
struct UnknownHeaderFlagError;

/// An Error that is thrown when a token parsing function
/// recognizes that the token vector index is out of bounds.
#[derive(Error, Debug)]
#[error("token index out of bounds")]
struct TokenIndexError;

/// An Error that is thrown when a token parsing function
/// recognizes that the token vector index is out of bounds.
#[derive(Error, Debug)]
#[error("unable to find a float-like substring")]
struct FloatRegexError;

/// The alphabet of the sequences represented in a P7HMM.
#[derive(Default)]
pub enum P7Alphabet {
    Amino,
    Dna,
    Rna,
    #[default]
    AlphabetNotSet,
}

/// Represents the header of the P7HMM.
#[derive(Default)]
pub struct Header {
    pub name: String,
    pub version: String,
    pub accession_number: String,
    pub description: String,
    pub date: String,
    pub command_line_history: String,
    pub has_reference_annotation: bool,
    pub has_model_mask: bool,
    pub has_consensus_residue: bool,
    pub has_consensus_structure: bool,
    pub has_map_annotation: bool,
    pub model_length: usize,
    pub max_length: usize,
    pub checksum: usize,
    pub num_sequences: usize,
    pub effective_num_sequences: f32,
    pub gathering_thresholds: [f32; 2],
    pub trusted_cutoffs: [f32; 2],
    pub noise_cutoffs: [f32; 2],
    pub alphabet: P7Alphabet,
}

/// This defines statistical scoring parameters for different pipeline stages.
#[derive(Default)]
pub struct Stats {
    pub msv_gumble_mu: f32,
    pub msv_gumble_lambda: f32,
    pub viterbi_gumble_mu: f32,
    pub viterbi_gumble_lambda: f32,
    pub forward_tau: f32,
    pub forward_lambda: f32,
}

/// The score model.
#[derive(Default)]
pub struct Model {
    pub composition: Vec<f32>,
    // TODO: optimize these as flat vectors
    pub match_probabilities: Vec<Vec<f32>>,
    pub insert_probabilities: Vec<Vec<f32>>,
    pub transition_probabilities: Vec<Vec<f32>>,
    // these remaining fields are the optional fields on the match emission line
    pub map_annotations: Vec<usize>,
    pub consensus_residues: String,
    pub reference_annotation: String,
    pub model_mask: Vec<bool>,
    pub consensus_structure: String,
}

/// The data that describes a complete P7HMM.
#[derive(Default)]
pub struct Hmm {
    pub header: Header,
    pub stats: Stats,
    pub model: Model,
}

impl Hmm {
    pub fn new() -> Self {
        let mut hmm = Hmm::default();

        hmm.model
            .match_probabilities
            .push(vec![0.0; MAX_DEGENERATE_ALPHABET_SIZE]);
        hmm.model
            .insert_probabilities
            .push(vec![0.0; MAX_DEGENERATE_ALPHABET_SIZE]);

        hmm.model.transition_probabilities.push(vec![0.0f32; 7]);

        hmm
    }
}

pub fn parse_hmms_from_p7hmm_file<R: AsRef<Path> + Display>(path: R) -> Result<Vec<Hmm>> {
    let phmm_file = File::open(&path)?;
    let mut phmm_lines = BufReader::new(phmm_file).lines();

    let mut line_number: usize = 0;
    let mut current_model_line_number: usize = 0;

    let mut parser_state = ParserState::Idle;
    let mut body_parser_state = ModelParserState::MatchEmissions;

    let mut hmm_list: Vec<Hmm> = vec![Hmm::new()];
    let mut current_hmm = hmm_list.get_mut(0).unwrap();
    let mut hmm_count: usize = 0;

    while let Some(Ok(line)) = phmm_lines.next() {
        let tokens: Vec<&str> = line.split_whitespace().collect();
        let flag: &str = get_token_as_str(&tokens, 0)?;

        line_number += 1;

        if line.trim().is_empty() {
            continue;
        }

        match parser_state {
            // if we're idle, we're searching for the next header
            ParserState::Idle => match flag {
                P7_HEADER_FORMAT_FLAG => {
                    if hmm_count > 0 {
                        hmm_list.push(Hmm::new());
                        current_hmm = hmm_list.last_mut().unwrap();
                    }

                    current_hmm.header.version =
                        get_joined_tokens(&tokens, 1).with_context(|| "")?;
                    parser_state = ParserState::Header;
                }
                _ => panic!("unknown flag: {}", flag),
            },
            ParserState::Header => {
                let error_context = || {
                    format!(
                        "failed to parse p7hmm file header: {}\n         on line: {}\n         with flag: {}",
                        &path.to_string(),
                        line_number,
                        flag
                    )
                };
                match flag {
                    P7_HEADER_NAME_FLAG => {
                        current_hmm.header.name =
                            get_token_as_string(&tokens, 1).with_context(error_context)?;
                    }
                    P7_HEADER_ACCESSION_FLAG => {
                        current_hmm.header.accession_number =
                            get_token_as_string(&tokens, 1).with_context(error_context)?;
                    }
                    P7_HEADER_DESCRIPTION_FLAG => {
                        current_hmm.header.description =
                            get_token_as_string(&tokens, 1).with_context(error_context)?;
                    }
                    P7_HEADER_LENGTH_FLAG => {
                        current_hmm.header.model_length =
                            get_token_as_usize(&tokens, 1).with_context(error_context)?;
                    }
                    P7_HEADER_MAXL_FLAG => {
                        current_hmm.header.max_length =
                            get_token_as_usize(&tokens, 1).with_context(error_context)?;
                    }
                    P7_HEADER_ALPHABET_FLAG => {
                        let value = get_token_as_str(&tokens, 1).with_context(error_context)?;
                        current_hmm.header.alphabet = match value {
                            "amino" => P7Alphabet::Amino,
                            "dna" => P7Alphabet::Dna,
                            "rna" => P7Alphabet::Rna,
                            _ => panic!("unknown alphabet"),
                        }
                    }
                    P7_HEADER_REFERENCE_FLAG => {
                        let value = get_token_as_str(&tokens, 1).with_context(error_context)?;
                        current_hmm.header.has_reference_annotation = match value {
                            "yes" => true,
                            "no" => false,
                            _ => panic!("unknown reference annotation value"),
                        }
                    }
                    P7_HEADER_MASK_FLAG => {
                        let value = get_token_as_str(&tokens, 1).with_context(error_context)?;
                        current_hmm.header.has_model_mask = match value {
                            "yes" => true,
                            "no" => false,
                            _ => panic!("unknown model mask value"),
                        }
                    }
                    P7_HEADER_CONSENSUS_RESIDUE_FLAG => {
                        let value = get_token_as_str(&tokens, 1).with_context(error_context)?;
                        current_hmm.header.has_consensus_residue = match value {
                            "yes" => true,
                            "no" => false,
                            _ => panic!("unknown consensus residue value"),
                        }
                    }
                    P7_HEADER_CONSENSUS_STRUCTURE_FLAG => {
                        let value = get_token_as_str(&tokens, 1).with_context(error_context)?;
                        current_hmm.header.has_consensus_structure = match value {
                            "yes" => true,
                            "no" => false,
                            _ => panic!("unknown consensus structure value"),
                        }
                    }
                    P7_HEADER_MAP_FLAG => {
                        let value = get_token_as_str(&tokens, 1).with_context(error_context)?;
                        current_hmm.header.has_map_annotation = match value {
                            "yes" => true,
                            "no" => false,
                            _ => panic!("unknown map annotation value"),
                        }
                    }
                    P7_HEADER_DATE_FLAG => {
                        current_hmm.header.date =
                            get_joined_tokens(&tokens, 1).with_context(error_context)?;
                    }
                    P7_HEADER_COMMAND_FLAG => {
                        // TODO: this is wrong I think, since the command line history might be multiple lines?
                        current_hmm.header.command_line_history =
                            get_joined_tokens(&tokens, 1).with_context(error_context)?;
                    }
                    P7_HEADER_NSEQ_FLAG => {
                        current_hmm.header.num_sequences =
                            get_token_as_usize(&tokens, 1).with_context(error_context)?;
                    }
                    P7_HEADER_EFFN_FLAG => {
                        current_hmm.header.effective_num_sequences =
                            get_token_as_f32(&tokens, 1).with_context(error_context)?;
                    }
                    P7_HEADER_CHECKSUM_FLAG => {
                        current_hmm.header.checksum =
                            get_token_as_usize(&tokens, 1).with_context(error_context)?;
                    }
                    P7_HEADER_GATHERING_FLAG => {
                        current_hmm.header.gathering_thresholds = [
                            get_token_as_f32(&tokens, 1).with_context(error_context)?,
                            get_token_as_f32(&tokens, 2).with_context(error_context)?,
                        ]
                    }
                    P7_HEADER_TRUSTED_FLAG => {
                        current_hmm.header.trusted_cutoffs = [
                            get_token_as_f32(&tokens, 1).with_context(error_context)?,
                            get_token_as_f32(&tokens, 2).with_context(error_context)?,
                        ]
                    }
                    P7_HEADER_NOISE_FLAG => {
                        current_hmm.header.noise_cutoffs = [
                            get_token_as_f32(&tokens, 1).with_context(error_context)?,
                            get_token_as_f32(&tokens, 2).with_context(error_context)?,
                        ]
                    }
                    P7_HEADER_STATS_FLAG => {
                        let mu_or_tau = get_token_as_f32(&tokens, 3).with_context(error_context)?;
                        let lambda = get_token_as_f32(&tokens, 4).with_context(error_context)?;
                        let sub_flag = get_token_as_str(&tokens, 1).with_context(error_context)?;
                        match sub_flag {
                            "LOCAL" => {
                                let sub_sub_flag =
                                    get_token_as_str(&tokens, 2).with_context(error_context)?;
                                match sub_sub_flag {
                                    P7_HEADER_STATS_MSV_FLAG => {
                                        current_hmm.stats.msv_gumble_mu = mu_or_tau;
                                        current_hmm.stats.msv_gumble_lambda = lambda;
                                    }
                                    P7_HEADER_STATS_VITERBI_FLAG => {
                                        current_hmm.stats.viterbi_gumble_mu = mu_or_tau;
                                        current_hmm.stats.viterbi_gumble_lambda = lambda;
                                    }
                                    P7_HEADER_STATS_FORWARD_FLAG => {
                                        current_hmm.stats.forward_tau = mu_or_tau;
                                        current_hmm.stats.forward_lambda = lambda;
                                    }
                                    _ => panic!("unknown header stats value (2)"),
                                }
                            }
                            _ => panic!("unknown stats value (1)"),
                        };
                    }
                    P7_BODY_HMM_MODEL_START_FLAG => {
                        // TODO: check against set alphabet?
                        parser_state = ParserState::ModelHead;
                    }
                    _ => {
                        // TODO: think about this?
                        // currently no-op for unknown flag
                        // println!("unknown flag: {}", flag)
                    }
                }
            }
            ParserState::ModelHead => {
                let error_context = || {
                    format!(
                        "failed to parse p7hmm file: {}\n         on line: {}\n",
                        &path.to_string(),
                        line_number
                    )
                };
                match body_parser_state {
                    ModelParserState::MatchEmissions => match flag {
                        P7_BODY_COMPO_FLAG => {
                            current_hmm.model.composition =
                                get_tokens_as_probability_vec(&tokens, 1, 21)
                                    .with_context(error_context)?;
                            body_parser_state = ModelParserState::InsertEmissions;
                        }
                        _ => {
                            // no-op for the the transition header
                        }
                    },
                    ModelParserState::InsertEmissions => {
                        current_hmm.model.insert_probabilities[0] =
                            get_tokens_as_probability_vec(&tokens, 0, 20)
                                .with_context(error_context)?;
                        body_parser_state = ModelParserState::StateTransitions;
                    }
                    ModelParserState::StateTransitions => {
                        current_hmm.model.transition_probabilities[0] =
                            get_tokens_as_probability_vec(&tokens, 0, 7)
                                .with_context(error_context)?;
                        current_model_line_number = 1;
                        parser_state = ParserState::ModelBody;
                        body_parser_state = ModelParserState::MatchEmissions;
                    }
                }
            }
            ParserState::ModelBody => {
                let error_context = || {
                    format!(
                        "failed to parse p7hmm file body: {}\n         on line: {}\n",
                        &path.to_string(),
                        line_number,
                    )
                };

                match flag {
                    P7_BODY_END_FLAG => {
                        // TODO: error handling to validate the parsed models
                        // for now let's at least just try to make sure we parsed the scores correctly
                        // we want it to be model_length + 1 since we are doing 1-indexing
                        assert_eq!(
                            current_hmm.model.match_probabilities.len(),
                            current_hmm.header.model_length + 1 as usize
                        );
                        assert_eq!(
                            current_hmm.model.insert_probabilities.len(),
                            current_hmm.header.model_length + 1 as usize
                        );

                        hmm_count += 1;
                        parser_state = ParserState::Idle;
                    }
                    _ => match body_parser_state {
                        ModelParserState::MatchEmissions => {
                            let line_number_flag = current_model_line_number.to_string();

                            if flag == line_number_flag {
                                current_hmm.model.match_probabilities.push(
                                    get_tokens_as_probability_vec(&tokens, 1, 21)
                                        .with_context(error_context)?,
                                );

                                body_parser_state = ModelParserState::InsertEmissions;
                            } else {
                                // TODO: error for mismatched line number for node
                                println!("{}", line);
                                panic!("mismatched line number");
                            }
                        }
                        ModelParserState::InsertEmissions => {
                            current_hmm.model.insert_probabilities.push(
                                get_tokens_as_probability_vec(&tokens, 0, 20)
                                    .with_context(error_context)?,
                            );
                            body_parser_state = ModelParserState::StateTransitions;
                        }
                        ModelParserState::StateTransitions => {
                            current_hmm.model.transition_probabilities.push(
                                get_tokens_as_probability_vec(&tokens, 0, 7)
                                    .with_context(error_context)?,
                            );
                            current_model_line_number += 1;
                            body_parser_state = ModelParserState::MatchEmissions;
                        }
                    },
                }
            }
        }
    }
    Ok(hmm_list)
}

fn token_index_check(tokens: &Vec<&str>, idx: usize) -> Result<()> {
    if tokens.len() + 1 < idx {
        return Err(TokenIndexError.into());
    }
    Ok(())
}

fn get_token_as_str<'a>(tokens: &'a Vec<&str>, idx: usize) -> Result<&'a str> {
    token_index_check(tokens, idx)?;
    Ok(tokens[idx])
}

fn get_token_as_string(tokens: &Vec<&str>, idx: usize) -> Result<String> {
    token_index_check(tokens, idx)?;
    Ok(String::from(tokens[idx]))
}

fn get_joined_tokens(tokens: &Vec<&str>, idx: usize) -> Result<String> {
    token_index_check(tokens, idx)?;
    Ok(tokens[idx..].join(" "))
}

fn get_token_as_f32(tokens: &Vec<&str>, idx: usize) -> Result<f32> {
    token_index_check(tokens, idx)?;

    // * means infinity in the spec
    if tokens[idx] == "*" {
        return Ok(0.0);
    }

    let float_str = match FLOAT_RE.find(tokens[idx]) {
        Some(str) => str,
        None => {
            return Err(FloatRegexError)
                .with_context(|| format!("failed to parse token \"{}\" as f32", tokens[idx]));
        }
    };

    float_str
        .as_str()
        .parse::<f32>()
        .with_context(|| format!("failed to parse token \"{}\" as f32", tokens[idx]))
}

fn get_token_as_usize(tokens: &Vec<&str>, idx: usize) -> Result<usize> {
    token_index_check(tokens, idx)?;
    tokens[idx]
        .parse::<usize>()
        .with_context(|| format!("failed to parse token \"{}\" as f32", tokens[idx]))
}

/// Get a string token as a probability, which equates to negating and exponentiating the float.
///
/// This is because in the HMMER3/f P7HMM spec, emissions and transitions are written as -ln(P).
fn get_token_as_probability(tokens: &Vec<&str>, idx: usize) -> Result<f32> {
    let float = get_token_as_f32(tokens, idx)?;
    Ok((-float).exp())
}

/// Get a collection of string tokens as a vector of probabilities, which equates
/// to negating and exponentiating the floats.
///
/// This is because in the HMMER3/f P7HMM spec, emissions and transitions are written as -ln(P).
fn get_tokens_as_probability_vec(tokens: &Vec<&str>, start: usize, end: usize) -> Result<Vec<f32>> {
    let mut vec: Vec<f32> = vec![];
    for i in start..end {
        vec.push(get_token_as_probability(&tokens, i)?)
    }
    Ok(vec)
}
