use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use crate::alphabet::{AMINO_ALPHABET, AMINO_BACKGROUND_FREQUENCIES};
use crate::structs::Profile;
use anyhow::{Context, Result};
use lazy_static::lazy_static;
use regex::Regex;
use thiserror::Error;

use self::constants::{
    HMM_DELETE_TO_DELETE, HMM_DELETE_TO_MATCH, HMM_INSERT_TO_INSERT, HMM_INSERT_TO_MATCH,
    HMM_MATCH_TO_DELETE, HMM_MATCH_TO_INSERT, HMM_MATCH_TO_MATCH,
};

use super::Sequence;

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
const P7_HEADER_FORMAT_FLAG: &str = "HMMER3/f";
const P7_HEADER_NAME_FLAG: &str = "NAME";
const P7_HEADER_ACCESSION_FLAG: &str = "ACC";
const P7_HEADER_DESCRIPTION_FLAG: &str = "DESC";
const P7_HEADER_LENGTH_FLAG: &str = "LENG";
const P7_HEADER_MAXL_FLAG: &str = "MAXL";
const P7_HEADER_ALPHABET_FLAG: &str = "ALPH";
const P7_HEADER_REFERENCE_FLAG: &str = "RF";
const P7_HEADER_MASK_FLAG: &str = "MM";
const P7_HEADER_CONSENSUS_RESIDUE_FLAG: &str = "CONS";
const P7_HEADER_CONSENSUS_STRUCTURE_FLAG: &str = "CS";
const P7_HEADER_MAP_FLAG: &str = "MAP";
const P7_HEADER_DATE_FLAG: &str = "DATE";
const P7_HEADER_COMMAND_FLAG: &str = "COM";
const P7_HEADER_NSEQ_FLAG: &str = "NSEQ";
const P7_HEADER_EFFN_FLAG: &str = "EFFN";
const P7_HEADER_CHECKSUM_FLAG: &str = "CKSUM";
const P7_HEADER_GATHERING_FLAG: &str = "GA";
const P7_HEADER_TRUSTED_FLAG: &str = "TC";
const P7_HEADER_NOISE_FLAG: &str = "NC";
const P7_HEADER_STATS_FLAG: &str = "STATS";
const P7_HEADER_STATS_MSV_FLAG: &str = "MSV";
const P7_HEADER_STATS_VITERBI_FLAG: &str = "VITERBI";
const P7_HEADER_STATS_FORWARD_FLAG: &str = "FORWARD";

const P7_BODY_HMM_MODEL_START_FLAG: &str = "HMM";
const P7_BODY_COMPO_FLAG: &str = "COMPO";
const P7_BODY_END_FLAG: &str = "//";

const BLOSUM_62_P_OPEN: f32 = 0.02;
const BLOSUM_62_P_EXTEND: f32 = 0.4;

const BLOSUM_62_CONDITIONAL_PROB: [[f32; 29]; 29] = [
    [
        0.2782, 0.0150, 0.0280, 0.0481, 0.0208, 0.0688, 0.0120, 0.0425, 0.0428, 0.0694, 0.0171,
        0.0217, 0.0348, 0.0285, 0.0390, 0.0929, 0.0535, 0.0666, 0.0044, 0.0159, 0.0000, 0.0498,
        0.1119, 0.0766, 0.0428, 0.0150, 1.0000, 0.0000, 0.0000,
    ],
    [
        0.0956, 0.3217, 0.0250, 0.0227, 0.0255, 0.0325, 0.0107, 0.0521, 0.0278, 0.0851, 0.0210,
        0.0194, 0.0226, 0.0185, 0.0253, 0.0603, 0.0477, 0.0595, 0.0073, 0.0195, 0.0000, 0.0444,
        0.1372, 0.0412, 0.0278, 0.3217, 1.0000, 0.0000, 0.0000,
    ],
    [
        0.0406, 0.0057, 0.3508, 0.1228, 0.0149, 0.0492, 0.0162, 0.0221, 0.0421, 0.0263, 0.0089,
        0.0554, 0.0342, 0.0385, 0.0279, 0.0665, 0.0383, 0.0252, 0.0031, 0.0114, 0.0000, 0.4062,
        0.0484, 0.1613, 0.0421, 0.0057, 1.0000, 0.0000, 0.0000,
    ],
    [
        0.0527, 0.0039, 0.0929, 0.3009, 0.0141, 0.0338, 0.0211, 0.0209, 0.0750, 0.0341, 0.0116,
        0.0381, 0.0323, 0.0686, 0.0497, 0.0628, 0.0361, 0.0327, 0.0040, 0.0148, 0.0000, 0.1309,
        0.0550, 0.3696, 0.0750, 0.0039, 1.0000, 0.0000, 0.0000,
    ],
    [
        0.0472, 0.0091, 0.0233, 0.0292, 0.3029, 0.0303, 0.0189, 0.0668, 0.0259, 0.1091, 0.0269,
        0.0181, 0.0153, 0.0173, 0.0236, 0.0410, 0.0324, 0.0555, 0.0178, 0.0894, 0.0000, 0.0414,
        0.1759, 0.0464, 0.0259, 0.0091, 1.0000, 0.0000, 0.0000,
    ],
    [
        0.0791, 0.0059, 0.0391, 0.0355, 0.0154, 0.4701, 0.0122, 0.0166, 0.0316, 0.0271, 0.0092,
        0.0416, 0.0257, 0.0210, 0.0288, 0.0686, 0.0287, 0.0260, 0.0061, 0.0118, 0.0000, 0.0807,
        0.0437, 0.0565, 0.0316, 0.0059, 1.0000, 0.0000, 0.0000,
    ],
    [
        0.0430, 0.0060, 0.0401, 0.0689, 0.0298, 0.0379, 0.3006, 0.0234, 0.0446, 0.0383, 0.0130,
        0.0587, 0.0263, 0.0408, 0.0558, 0.0512, 0.0295, 0.0267, 0.0062, 0.0592, 0.0000, 0.0988,
        0.0617, 0.1096, 0.0446, 0.0060, 1.0000, 0.0000, 0.0000,
    ],
    [
        0.0586, 0.0113, 0.0211, 0.0263, 0.0406, 0.0199, 0.0090, 0.2152, 0.0234, 0.1861, 0.0334,
        0.0163, 0.0190, 0.0156, 0.0213, 0.0370, 0.0402, 0.1787, 0.0045, 0.0226, 0.0000, 0.0374,
        0.4012, 0.0419, 0.0234, 0.0113, 1.0000, 0.0000, 0.0000,
    ],
    [
        0.0550, 0.0056, 0.0374, 0.0881, 0.0147, 0.0353, 0.0160, 0.0218, 0.2796, 0.0490, 0.0166,
        0.0398, 0.0337, 0.0522, 0.0980, 0.0656, 0.0378, 0.0342, 0.0042, 0.0155, 0.0000, 0.0771,
        0.0708, 0.1403, 0.2796, 0.0056, 1.0000, 0.0000, 0.0000,
    ],
    [
        0.0580, 0.0112, 0.0152, 0.0260, 0.0401, 0.0197, 0.0089, 0.1127, 0.0318, 0.3476, 0.0454,
        0.0161, 0.0188, 0.0212, 0.0290, 0.0366, 0.0398, 0.0936, 0.0061, 0.0224, 0.0000, 0.0313,
        0.4602, 0.0472, 0.0318, 0.0112, 1.0000, 0.0000, 0.0000,
    ],
    [
        0.0595, 0.0114, 0.0214, 0.0367, 0.0412, 0.0278, 0.0126, 0.0841, 0.0448, 0.1887, 0.1208,
        0.0228, 0.0265, 0.0410, 0.0408, 0.0516, 0.0408, 0.0960, 0.0086, 0.0229, 0.0000, 0.0441,
        0.2728, 0.0777, 0.0448, 0.0114, 1.0000, 0.0000, 0.0000,
    ],
    [
        0.0406, 0.0057, 0.0716, 0.0650, 0.0149, 0.0676, 0.0306, 0.0221, 0.0578, 0.0361, 0.0122,
        0.2716, 0.0249, 0.0385, 0.0526, 0.0914, 0.0526, 0.0252, 0.0031, 0.0157, 0.0000, 0.3432,
        0.0582, 0.1035, 0.0578, 0.0057, 1.0000, 0.0000, 0.0000,
    ],
    [
        0.0578, 0.0059, 0.0392, 0.0490, 0.0112, 0.0371, 0.0122, 0.0229, 0.0436, 0.0374, 0.0127,
        0.0221, 0.4505, 0.0290, 0.0289, 0.0501, 0.0396, 0.0359, 0.0032, 0.0118, 0.0000, 0.0613,
        0.0603, 0.0780, 0.0436, 0.0059, 1.0000, 0.0000, 0.0000,
    ],
    [
        0.0566, 0.0058, 0.0528, 0.1245, 0.0151, 0.0363, 0.0226, 0.0224, 0.0806, 0.0503, 0.0235,
        0.0409, 0.0347, 0.1914, 0.0734, 0.0674, 0.0388, 0.0352, 0.0060, 0.0218, 0.0000, 0.0937,
        0.0728, 0.3159, 0.0806, 0.0058, 1.0000, 0.0000, 0.0000,
    ],
    [
        0.0598, 0.0061, 0.0296, 0.0697, 0.0160, 0.0384, 0.0239, 0.0237, 0.1171, 0.0532, 0.0180,
        0.0432, 0.0267, 0.0567, 0.2766, 0.0519, 0.0410, 0.0271, 0.0046, 0.0168, 0.0000, 0.0728,
        0.0769, 0.1264, 0.1171, 0.0061, 1.0000, 0.0000, 0.0000,
    ],
    [
        0.1030, 0.0105, 0.0509, 0.0636, 0.0200, 0.0661, 0.0159, 0.0297, 0.0566, 0.0485, 0.0165,
        0.0542, 0.0334, 0.0376, 0.0374, 0.2319, 0.0707, 0.0339, 0.0042, 0.0153, 0.0000, 0.1051,
        0.0783, 0.1012, 0.0566, 0.0105, 1.0000, 0.0000, 0.0000,
    ],
    [
        0.0777, 0.0109, 0.0384, 0.0479, 0.0207, 0.0363, 0.0120, 0.0423, 0.0426, 0.0691, 0.0171,
        0.0409, 0.0346, 0.0284, 0.0388, 0.0926, 0.2614, 0.0664, 0.0060, 0.0159, 0.0000, 0.0793,
        0.1115, 0.0763, 0.0426, 0.0109, 1.0000, 0.0000, 0.0000,
    ],
    [
        0.0788, 0.0110, 0.0206, 0.0354, 0.0289, 0.0268, 0.0088, 0.1531, 0.0315, 0.1324, 0.0327,
        0.0160, 0.0256, 0.0209, 0.0208, 0.0362, 0.0540, 0.2401, 0.0044, 0.0221, 0.0000, 0.0366,
        0.2855, 0.0563, 0.0315, 0.0110, 1.0000, 0.0000, 0.0000,
    ],
    [
        0.0345, 0.0091, 0.0171, 0.0293, 0.0621, 0.0419, 0.0138, 0.0259, 0.0261, 0.0581, 0.0197,
        0.0132, 0.0154, 0.0238, 0.0237, 0.0300, 0.0326, 0.0295, 0.4290, 0.0654, 0.0000, 0.0303,
        0.0839, 0.0531, 0.0261, 0.0091, 1.0000, 0.0000, 0.0000,
    ],
    [
        0.0437, 0.0084, 0.0216, 0.0370, 0.1078, 0.0280, 0.0453, 0.0449, 0.0329, 0.0734, 0.0181,
        0.0230, 0.0195, 0.0301, 0.0300, 0.0379, 0.0300, 0.0513, 0.0226, 0.2947, 0.0000, 0.0445,
        0.1183, 0.0671, 0.0329, 0.0084, 1.0000, 0.0000, 0.0000,
    ],
    [
        0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000,
        0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000,
        0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000,
    ],
    [
        0.0406, 0.0057, 0.2290, 0.0976, 0.0149, 0.0572, 0.0225, 0.0221, 0.0490, 0.0306, 0.0104,
        0.1498, 0.0301, 0.0385, 0.0387, 0.0773, 0.0445, 0.0252, 0.0031, 0.0133, 0.0000, 0.3787,
        0.0527, 0.1361, 0.0490, 0.0057, 1.0000, 0.0000, 0.0000,
    ],
    [
        0.0582, 0.0112, 0.0174, 0.0261, 0.0403, 0.0198, 0.0090, 0.1513, 0.0286, 0.2866, 0.0409,
        0.0162, 0.0189, 0.0191, 0.0261, 0.0367, 0.0399, 0.1257, 0.0055, 0.0225, 0.0000, 0.0336,
        0.4380, 0.0452, 0.0286, 0.0112, 1.0000, 0.0000, 0.0000,
    ],
    [
        0.0541, 0.0046, 0.0786, 0.2383, 0.0144, 0.0347, 0.0216, 0.0214, 0.0770, 0.0399, 0.0158,
        0.0391, 0.0331, 0.1122, 0.0581, 0.0644, 0.0371, 0.0336, 0.0047, 0.0173, 0.0000, 0.1177,
        0.0613, 0.3505, 0.0770, 0.0046, 1.0000, 0.0000, 0.0000,
    ],
    [
        0.0550, 0.0056, 0.0374, 0.0881, 0.0147, 0.0353, 0.0160, 0.0218, 0.2796, 0.0490, 0.0166,
        0.0398, 0.0337, 0.0522, 0.0980, 0.0656, 0.0378, 0.0342, 0.0042, 0.0155, 0.0000, 0.0771,
        0.0708, 0.1403, 0.2796, 0.0056, 1.0000, 0.0000, 0.0000,
    ],
    [
        0.0956, 0.3217, 0.0250, 0.0227, 0.0255, 0.0325, 0.0107, 0.0521, 0.0278, 0.0851, 0.0210,
        0.0194, 0.0226, 0.0185, 0.0253, 0.0603, 0.0477, 0.0595, 0.0073, 0.0195, 0.0000, 0.0444,
        0.1372, 0.0412, 0.0278, 0.3217, 1.0000, 0.0000, 0.0000,
    ],
    [
        0.0796, 0.0125, 0.0550, 0.0728, 0.0351, 0.0692, 0.0222, 0.0577, 0.0619, 0.0953, 0.0229,
        0.0426, 0.0479, 0.0401, 0.0519, 0.0718, 0.0548, 0.0674, 0.0100, 0.0291, 0.0000, 0.0976,
        0.1531, 0.1129, 0.0619, 0.0125, 1.0000, 0.0000, 0.0000,
    ],
    [
        0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000,
        0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000,
        0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000,
    ],
    [
        0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000,
        0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000,
        0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000,
    ],
];

// this static regex is used to find float strings
lazy_static! {
    static ref FLOAT_RE: Regex = Regex::new(r"-*\d\.*\d*").unwrap();
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
#[derive(Default, Clone)]
pub enum Alphabet {
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
    pub alphabet: Alphabet,
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
            .push(vec![0.0; Profile::MAX_DEGENERATE_ALPHABET_SIZE]);
        hmm.model
            .insert_probabilities
            .push(vec![0.0; Profile::MAX_DEGENERATE_ALPHABET_SIZE]);

        hmm.model.transition_probabilities.push(vec![0.0f32; 7]);

        hmm
    }

    pub fn from_blosum_62_and_sequence(seq: &Sequence) -> Result<Self> {
        let mut match_probabilities = vec![vec![0.0; 20]; seq.length + 1];

        seq.digital_bytes
            .iter()
            .enumerate()
            .skip(1)
            .for_each(|(pos, &residue)| {
                match_probabilities[pos] =
                    BLOSUM_62_CONDITIONAL_PROB[residue as usize][0..AMINO_ALPHABET.len()].to_vec();
            });

        // the transition probabilities are uniform
        let mut transitions = vec![0.0; 8];
        transitions[HMM_MATCH_TO_MATCH] = 1.0 - 2.0 * BLOSUM_62_P_OPEN;
        transitions[HMM_MATCH_TO_INSERT] = BLOSUM_62_P_OPEN;
        transitions[HMM_MATCH_TO_DELETE] = BLOSUM_62_P_OPEN;
        transitions[HMM_INSERT_TO_MATCH] = 1.0 - BLOSUM_62_P_EXTEND;
        transitions[HMM_INSERT_TO_INSERT] = BLOSUM_62_P_EXTEND;
        transitions[HMM_DELETE_TO_MATCH] = 1.0 - BLOSUM_62_P_EXTEND;
        transitions[HMM_DELETE_TO_DELETE] = BLOSUM_62_P_EXTEND;

        let mut transition_probabilities = vec![transitions; seq.length + 1];
        // except for a few modifications to the last position
        transition_probabilities[seq.length][HMM_MATCH_TO_MATCH] = 1.0 - BLOSUM_62_P_OPEN;
        transition_probabilities[seq.length][HMM_MATCH_TO_DELETE] = 0.0;
        transition_probabilities[seq.length][HMM_DELETE_TO_MATCH] = 1.0;
        transition_probabilities[seq.length][HMM_DELETE_TO_DELETE] = 0.0;

        let insert_probabilities = vec![AMINO_BACKGROUND_FREQUENCIES.to_vec(); seq.length + 1];

        // OCCUPANCY

        let mut match_occupancy = vec![0.0f32; seq.length + 1];
        match_occupancy[1] = transition_probabilities[0][HMM_MATCH_TO_INSERT]
            + transition_probabilities[0][HMM_MATCH_TO_MATCH];

        (2..=seq.length)
            .map(|pos| (pos - 1, pos))
            .for_each(|(prev_pos, pos)| {
                let a = match_occupancy[prev_pos]
                    * (transition_probabilities[prev_pos][HMM_MATCH_TO_MATCH]
                        + transition_probabilities[prev_pos][HMM_MATCH_TO_INSERT]);

                let b = (1.0 - match_occupancy[prev_pos])
                    * transition_probabilities[prev_pos][HMM_DELETE_TO_MATCH];

                match_occupancy[pos] = a + b;
            });

        let mut insert_occupancy = vec![0.0f32; seq.length + 1];
        insert_occupancy[0] = transition_probabilities[0][HMM_MATCH_TO_INSERT]
            / transition_probabilities[0][HMM_INSERT_TO_MATCH];

        (1..=seq.length).for_each(|pos| {
            insert_occupancy[pos] = match_occupancy[pos]
                * transition_probabilities[pos][HMM_MATCH_TO_INSERT]
                / transition_probabilities[pos][HMM_INSERT_TO_MATCH]
        });

        // COMPOSITION

        let mut composition = vec![0.0f32; AMINO_ALPHABET.len()];

        let mut add_scaled_vec_fn = |probs, scalar| {
            composition
                .iter_mut()
                .zip(probs)
                .for_each(|(c, p)| *c += p * scalar)
        };

        add_scaled_vec_fn(&insert_probabilities[0], insert_occupancy[0]);
        (1..=seq.length).for_each(|pos| {
            add_scaled_vec_fn(&match_probabilities[pos], match_occupancy[pos]);
            add_scaled_vec_fn(&insert_probabilities[pos], insert_occupancy[pos]);
        });

        let composition_sum = composition.iter().sum::<f32>();
        composition.iter_mut().for_each(|c| *c /= composition_sum);

        // LAMBDA

        let total_entropy = match_probabilities
            .iter()
            .skip(1)
            .flat_map(|probs| probs.iter().zip(AMINO_BACKGROUND_FREQUENCIES))
            .map(|(p, f)| p * (p / f).log2())
            .sum::<f32>();

        let lambda = std::f32::consts::LN_2 + 1.44 / total_entropy;

        Ok(Self {
            header: Header {
                name: seq.name.clone(),
                model_length: seq.length,
                alphabet: Alphabet::Amino,
                ..Default::default()
            },
            stats: Stats {
                forward_lambda: lambda,
                msv_gumble_lambda: lambda,
                viterbi_gumble_lambda: lambda,
                ..Default::default()
            },
            model: Model {
                composition,
                match_probabilities,
                insert_probabilities,
                transition_probabilities,
                consensus_residues: String::from_utf8(seq.utf8_bytes[1..].to_vec())?,
                ..Default::default()
            },
        })
    }
}

// TODO: remove the Display trait bound
pub fn parse_hmms_from_p7hmm_file<P: AsRef<Path>>(path: P) -> Result<Vec<Hmm>> {
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
                        &path.as_ref().to_string_lossy(),
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
                            "amino" => Alphabet::Amino,
                            "dna" => Alphabet::Dna,
                            "rna" => Alphabet::Rna,
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
                        &path.as_ref().to_string_lossy(),
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
                        &path.as_ref().to_string_lossy(),
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
                            current_hmm.header.model_length + 1
                        );
                        assert_eq!(
                            current_hmm.model.insert_probabilities.len(),
                            current_hmm.header.model_length + 1
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
        vec.push(get_token_as_probability(tokens, i)?)
    }
    Ok(vec)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::structs::Sequence;

    #[test]
    fn test_hmm_from_blosum62() -> anyhow::Result<()> {
        let mut seq = Sequence::from_utf8(
            concat!(
                "GNLLVILVILRNKKLRTPTNIFLLNLAVADLLVLLLVLPFSLVYALLEGDWVFGEVLCKL",
                "VTALDVVNLTASILLLTAISIDRYLAIVKPLKYKRIRTKRRALVLILVVWVLALLLSLPP",
                "LLFSGTKTESAEKEETVCLIDFPEEESTWEVSYTLLLSVLGFLLPLLVILVCYVRILRTL",
                "RKSAKKEKSRKKKSARKERKALKTLLVVVVVFVLCWLPYFILLLLDSLLKECESEKLVET",
                "ALLITLLLAYVNSCLNPIIY"
            )
            .as_bytes(),
        )?;
        seq.name = ">7tm_1-consensus".to_string();

        let hmm = Hmm::from_blosum_62_and_sequence(&seq)?;

        let composition_corrrect = [
            0.075741, 0.016550, 0.039336, 0.060816, 0.038081, 0.043406, 0.014495, 0.077001,
            0.061813, 0.144704, 0.026453, 0.033726, 0.040502, 0.031336, 0.048901, 0.063073,
            0.056054, 0.085715, 0.012350, 0.029948,
        ];

        hmm.model
            .composition
            .iter()
            .zip(composition_corrrect)
            .for_each(|(a, b)| assert!((a - b).abs() < 1e-4));

        let lambda_correct = 0.703323;
        assert!((lambda_correct - hmm.stats.forward_lambda) < 1e-4,);

        Profile::new(&hmm);
        Ok(())
    }
}
