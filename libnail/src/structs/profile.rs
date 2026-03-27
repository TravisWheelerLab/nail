use std::cmp::Ordering;
use std::fmt::{self, Debug};
use std::fmt::{Display, Formatter};
use std::ops::Index;

use super::Sequence;
use crate::alphabet::{AminoUtilsDigital, UTF8_SPACE};
use crate::util::{LogAbuse, VecMath};
use crate::{
    align::{
        forward, null_one_score,
        structs::{BackgroundState, CoreState, DpMatrixSparse, RowBounds, Trace},
    },
    alphabet::{Alphabet, AminoAcid, AMINO_ALPHABET_WITH_DEGENERATE, AMINO_BACKGROUND_FREQUENCIES},
    util::mean_relative_entropy,
};

use anyhow::{anyhow, Context};
use datasize::DataSize;
use rand::SeedableRng;
use rand_pcg::Pcg64;

impl AsRef<Profile> for Profile {
    fn as_ref(&self) -> &Profile {
        self
    }
}

mod blosum62 {
    pub const P_OPEN: f32 = 0.02;
    pub const P_EXTEND: f32 = 0.4;

    pub const CONDITIONAL_PROB: [[f32; 29]; 29] = [
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
}

pub struct Emission(pub CoreState, pub AminoAcid);
pub struct CoreToCore(pub CoreState, pub CoreState);

impl Debug for CoreToCore {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        write!(f, "{:?} -> {:?}", self.0, self.1)
    }
}
pub struct BackgroundLoop(pub BackgroundState);
pub struct CoreEntry(pub CoreState);

/// The Cantor pairing function:
///   π(a, b) = (a + b) * (a + b + 1) / 2 + b
///
/// This is used to map each pair of core states
/// (i.e. a transition) to a unique integer.
//
// this happens to work out (almost) perfectly for mapping
// pairs of states for valid core transitions:
//   +----------+--------+--------------------------------+---------+
//   | Symbolic |  Pair  |          Computation           | π(a, b) |
//   +----------+--------+--------------------------------+---------+
//   |  (M, M)  | (0, 0) | (0 + 0) * (0 + 0 + 1) / 2 + 0  |    0    |
//   |  (I, M)  | (1, 0) | (1 + 0) * (1 + 0 + 1) / 2 + 0  |    1    |
//   |  (M, I)  | (0, 1) | (0 + 1) * (0 + 1 + 1) / 2 + 1  |    2    |
//   |  (D, M)  | (2, 0) | (2 + 0) * (2 + 0 + 1) / 2 + 0  |    3    |
//   |  (I, I)  | (1, 1) | (1 + 1) * (1 + 1 + 1) / 2 + 1  |    4    |
//   |  (M, D)  | (0, 2) | (0 + 2) * (0 + 2 + 1) / 2 + 2  |    5    |
//   |  (D, D)  | (2, 2) | (2 + 2) * (2 + 2 + 1) / 2 + 2  |   12    |
//   |  (D, I)  | (2, 1) | (2 + 1) * (2 + 1 + 1) / 2 + 1  |    7    |
//   |  (I, D)  | (1, 2) | (1 + 2) * (1 + 2 + 1) / 2 + 2  |    8    |
//   +----------+--------+--------------------------------+---------+
fn cantor(a: usize, b: usize) -> usize {
    // debug guard rails to catch D->I or I->D
    debug_assert!(a < 3);
    debug_assert!(b < 3);
    (a + b) * (a + b + 1) / 2 + b
}

mod serialize {
    use std::io::{Read, Write};

    use super::Profile;

    #[repr(u8)]
    enum Version {
        V1 = 0u8,
    }

    fn read_u32_le<R: Read>(r: &mut R) -> anyhow::Result<u32> {
        let mut buf = [0u8; 4];
        r.read_exact(&mut buf)?;
        Ok(u32::from_le_bytes(buf))
    }

    fn read_f32_le<R: Read>(r: &mut R) -> anyhow::Result<f32> {
        let mut buf = [0u8; 4];
        r.read_exact(&mut buf)?;
        Ok(f32::from_le_bytes(buf))
    }

    impl Profile {
        pub fn serialize<W: Write>(&self, out: &mut W) -> std::io::Result<()> {
            // [HEADER] | 2 bytes <custom>
            out.write_all(&[Version::V1 as u8, 0])?;

            // [ALPHABET] | 1 byte <u8>
            out.write_all(&[self.alphabet as u8])?;

            // [MODEL LEN] | 4 bytes <u32>
            out.write_all(&(self.length as u32).to_le_bytes())?;

            // [NAME LEN] | 4 bytes <u32>
            out.write_all(&(self.name.len() as u32).to_le_bytes())?;

            // [ACCESSION LEN] | 4 bytes <u32>
            out.write_all(&(self.accession.len() as u32).to_le_bytes())?;

            // [TAU] | 4 bytes <f32>
            out.write_all(&self.fwd_tau.to_le_bytes())?;

            // [LAMBDA] | 4 bytes <f32>
            out.write_all(&self.fwd_lambda.to_le_bytes())?;

            // [NAME] | N bytes <u8>(UTF8)
            out.write_all(self.name.as_bytes())?;

            // [ACCESSION] | A Bytes <u8>(UTF8)
            out.write_all(self.accession.as_bytes())?;

            // [CONSENSUS] | M bytes <u8>(UTF8)
            out.write_all(&self.consensus_seq_bytes_utf8)?;

            // [TRANSITIONS] | [ 4 * 8 ] * ( M + 1 ) bytes
            for v in &self.core_transitions {
                for f in v {
                    out.write_all(&f.to_le_bytes())?;
                }
            }

            // [EMISSIONS] | [ 4 * 29 ] * ( M + 2 ) bytes
            for v in &self.emission_scores[0] {
                for f in v {
                    out.write_all(&f.to_le_bytes())?;
                }
            }

            Ok(())
        }

        pub fn deserialize<R: Read>(mut buf: R) -> anyhow::Result<Self> {
            let mut header = [0u8; 2];
            buf.read_exact(&mut header)?;

            let mut alphabet = [0u8; 1];
            buf.read_exact(&mut alphabet)?;

            let length = read_u32_le(&mut buf)? as usize;
            let name_len = read_u32_le(&mut buf)? as usize;
            let accession_len = read_u32_le(&mut buf)? as usize;

            let fwd_tau = read_f32_le(&mut buf)?;
            let fwd_lambda = read_f32_le(&mut buf)?;

            let mut name = vec![0u8; name_len];
            buf.read_exact(&mut name)?;

            let mut accession = vec![0u8; accession_len];
            buf.read_exact(&mut accession)?;

            let mut consensus = vec![0u8; length + 1];
            buf.read_exact(&mut consensus)?;

            let mut core_transitions = Vec::with_capacity(length + 1);
            for _ in 0..(length + 1) {
                let mut row = [0f32; 8];
                for f in &mut row {
                    *f = read_f32_le(&mut buf)?;
                }
                core_transitions.push(row);
            }

            let mut mat_emissions = Vec::with_capacity(length + 2);
            for _ in 0..(length + 2) {
                let mut row = [0f32; 29];
                for f in &mut row {
                    *f = read_f32_le(&mut buf)?;
                }
                mat_emissions.push(row);
            }

            Ok(Self {
                alphabet: alphabet[0].into(),
                length,
                name: String::from_utf8(name).unwrap(),
                accession: String::from_utf8(accession).unwrap(),
                consensus_seq_bytes_utf8: consensus,
                fwd_tau,
                fwd_lambda,
                core_transitions,
                emission_scores: [mat_emissions, vec![]],
                ..Default::default()
            })
        }
    }
}

#[derive(Clone, Default, DataSize)]
pub struct Profile {
    /// The name of the profile
    pub name: String,
    /// The accession number of the profile
    pub accession: String,
    /// Model length (number of nodes)
    pub length: usize,
    /// Current target sequence length
    pub target_length: usize,
    /// Calculated upper bound on max sequence length
    pub max_length: usize,
    /// Core model transition scores
    pub core_transitions: Vec<[f32; 8]>,
    /// Core model entry transition scores
    pub entry_transitions: Vec<f32>,
    /// Match and insert emission scores
    pub emission_scores: [Vec<[f32; Profile::MAX_DEGENERATE_ALPHABET_SIZE]>; 2],
    /// Transitions from special states (E, N, B, J, C)
    pub special_transitions: [[f32; 2]; 5],
    /// The expected number of times that the J state is used
    pub expected_j_uses: f32,
    /// The profile's consensus sequence
    pub consensus_seq_bytes_utf8: Vec<u8>,
    /// The sequence alphabet
    pub alphabet: Alphabet,
    pub fwd_tau: f32,
    pub fwd_lambda: f32,
}

impl Display for Profile {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        writeln!(f, "NAME  {}", self.name)?;
        writeln!(f, "LENG  {}", self.length)?;
        writeln!(f, "ALPH  {}", self.alphabet)?;
        writeln!(f, "TAU   {:.6}", self.fwd_tau)?;
        writeln!(f, "LAMB  {:.6}", self.fwd_lambda)?;

        let pw = self.length.to_string().len();

        write!(f, "{:W$}", "", W = pw + 2)?;
        self.alphabet
            .canonical_iter()
            .try_for_each(|c| write!(f, "  {:>6}", c.to_utf8_byte_amino() as char))?;
        writeln!(f)?;

        write!(f, "{:W$}", "", W = pw + 2)?;
        writeln!(
            f,
            "  {:>6}  {:>6}  {:>6}  {:>6}  {:>6}  {:>6}  {:>6}  {:>6}",
            "m->m", "i->m", "m->i", "d->m", "i->i", "m->d", "d->d", "b->m",
        )?;

        for i in 0..=self.length {
            write!(
                f,
                "{} {i:>W$}",
                self.consensus_seq_bytes_utf8[i] as char,
                W = pw
            )?;

            for s in &self.emission_scores[0][i] {
                write!(f, "  {s:6.3}")?;
            }
            writeln!(f)?;

            write!(f, "{:W$}", "", W = pw + 2)?;
            for s in &self.core_transitions[i] {
                write!(f, "  {s:6.3}")?;
            }
            writeln!(f)?;
            writeln!(f)?;
        }

        Ok(())
    }
}

impl Index<CoreToCore> for Profile {
    type Output = f32;

    fn index(&self, index: CoreToCore) -> &Self::Output {
        let profile_idx = index.0.profile_idx();
        let transition_idx = cantor(index.0.state_idx(), index.1.state_idx()).min(6);
        &self.core_transitions[profile_idx][transition_idx]
    }
}

impl Index<BackgroundLoop> for Profile {
    type Output = f32;

    fn index(&self, index: BackgroundLoop) -> &Self::Output {
        &self.special_transitions[index.0 as usize][Self::SPECIAL_LOOP_IDX]
    }
}

impl Index<CoreEntry> for Profile {
    type Output = f32;

    fn index(&self, index: CoreEntry) -> &Self::Output {
        let profile_idx = index.0.profile_idx();
        &self.entry_transitions[profile_idx]
    }
}

impl Index<Emission> for Profile {
    type Output = f32;

    fn index(&self, index: Emission) -> &Self::Output {
        let state_idx = index.0.state_idx();
        let profile_idx = index.0.profile_idx();
        let residue_idx = index.1 as usize;
        &self.emission_scores[state_idx][profile_idx][residue_idx]
    }
}

#[repr(usize)]
pub enum Transition {
    MM = 0,
    IM = 1,
    MI = 2,
    DM = 3,
    II = 4,
    MD = 5,
    DD = 6,
    BM = 7,
}

#[derive(Default)]
pub struct ProfileBuilder {
    name: Option<String>,
    accession: Option<String>,
    length: Option<usize>,
    mat_emissions: Vec<Option<[f32; Profile::MAX_ALPHABET_SIZE]>>,
    transitions: Vec<Option<[f32; 7]>>,
    fwd_tau: Option<f32>,
    fwd_lambda: Option<f32>,
}

impl ProfileBuilder {
    pub fn name(&mut self, name: String) -> &mut Self {
        self.name = Some(name);
        self
    }

    pub fn accession(&mut self, accession: String) -> &mut Self {
        self.accession = Some(accession);
        self
    }

    pub fn length(&mut self, length: usize) -> &mut Self {
        self.length = Some(length);
        self.mat_emissions.resize(length + 1, None);
        self.transitions.resize(length + 1, None);
        self
    }

    pub fn mat_emission(&mut self, pos: usize, probs: [f32; 20]) -> &mut Self {
        self.mat_emissions[pos] = Some(probs);
        self
    }

    pub fn transition(&mut self, pos: usize, probs: [f32; 7]) -> &mut Self {
        self.transitions[pos] = Some(probs);
        self
    }

    pub fn fwd_tau(&mut self, fwd_tau: f32) -> &mut Self {
        self.fwd_tau = Some(fwd_tau);
        self
    }

    pub fn fwd_lambda(&mut self, fwd_lambda: f32) -> &mut Self {
        self.fwd_lambda = Some(fwd_lambda);
        self
    }

    pub fn build(mut self) -> anyhow::Result<Profile> {
        let name = self.name.ok_or(anyhow!("missing: name"))?;
        let length = self.length.ok_or(anyhow!("missing: length"))?;
        let accession = self.accession.unwrap_or_default();
        let fwd_tau = self.fwd_tau.ok_or(anyhow!("missing: tau"))?;
        let fwd_lambda = self.fwd_lambda.ok_or(anyhow!("missing: tau"))?;

        self.mat_emissions[0] = Some([-f32::INFINITY; 20]);

        let mat_emissions = self
            .mat_emissions
            .into_iter()
            .enumerate()
            .map(|(i, v)| v.ok_or(anyhow!("missing: match emissions {}", i)))
            .collect::<anyhow::Result<Vec<_>>>()?;

        let transitions = self
            .transitions
            .into_iter()
            .enumerate()
            .map(|(i, v)| v.ok_or(anyhow!("missing: transitions {}", i + 1)))
            .collect::<anyhow::Result<Vec<_>>>()?;

        let mut prf = Profile {
            name,
            accession,
            length,
            target_length: 0,
            max_length: 0,
            core_transitions: vec![[0.0; Profile::NUM_STATE_TRANSITIONS]; length + 1],
            entry_transitions: vec![0.0; length + 1],
            emission_scores: [
                // +1 for left-pad & +1 for right-pad
                vec![[-f32::INFINITY; Profile::MAX_DEGENERATE_ALPHABET_SIZE]; length + 2],
                vec![[-f32::INFINITY; Profile::MAX_DEGENERATE_ALPHABET_SIZE]; length + 2],
            ],
            special_transitions: [[0.0; 2]; 5],
            expected_j_uses: 0.0,
            consensus_seq_bytes_utf8: vec![UTF8_SPACE; length + 1],
            alphabet: Alphabet::Amino,
            fwd_tau,
            fwd_lambda,
        };

        for state in 0..Profile::NUM_STATE_TRANSITIONS {
            prf.core_transitions[0][state] = -f32::INFINITY;
        }

        // occupancy
        let mut occupancy = vec![0.0; prf.length + 1];
        occupancy[1] =
            transitions[0][Transition::MI as usize] + transitions[0][Transition::MM as usize];
        for profile_idx in 2..=prf.length {
            // the occupancy of a model position is the
            // sum of the following two probabilities:
            occupancy[profile_idx] = (
                // the occupancy probability of the previous position
                occupancy[profile_idx - 1]
                    // multiplied by the sum of the transitions to "occupying" states
                    * (transitions[profile_idx - 1][Transition::MM as usize]
                        + transitions[profile_idx - 1][Transition::MI as usize])
            ) + (
                // the complement of the occupancy of the previous position
                (1.0 - occupancy[profile_idx - 1])
                    // multiplied by the transition to a match state
                    //   ** since there's no delete to insert transition **
                    * transitions[profile_idx - 1][Transition::DM as usize]
            )
        }

        let occupancy_sum: f32 = (1..=prf.length).fold(0.0, |acc, profile_idx| {
            acc + occupancy[profile_idx] * (prf.length - profile_idx + 1) as f32
        });

        // the model entry distribution is essentially the normalized occupancy
        (1..=prf.length).for_each(|profile_idx| {
            prf.core_transitions[profile_idx - 1][Transition::BM as usize] =
                (occupancy[profile_idx] / occupancy_sum).ln();

            prf.entry_transitions[profile_idx - 1] = (occupancy[profile_idx] / occupancy_sum).ln();
        });

        // these settings are for the non-multi-hit mode
        // N, C, and J transitions are set later by length config
        prf.special_transitions[Profile::E_IDX][Profile::SPECIAL_MOVE_IDX] = 0.0;
        prf.special_transitions[Profile::E_IDX][Profile::SPECIAL_LOOP_IDX] = -f32::INFINITY;
        prf.expected_j_uses = 0.0;

        // transition scores
        prf.core_transitions
            .iter_mut()
            .zip(transitions)
            .skip(1)
            .for_each(|(scores, probs)| {
                scores
                    .iter_mut()
                    .zip(probs)
                    .for_each(|(s, p)| *s = p.ln_or_inf())
            });

        prf.consensus_seq_bytes_utf8
            .iter_mut()
            .zip(mat_emissions.iter())
            .skip(1)
            .for_each(|(c, probs)| {
                // choose the residue with the highest
                // match emission probability
                //
                // TODO: this should not be the case
                //       for single sequence models
                //       (it should just be the seq)
                //
                // note: unwrap() because we check to
                //       make sure these exist earlier
                let residue = probs.argmax().unwrap();

                if probs[residue] > 0.5 {
                    *c = residue.to_utf8_byte_amino()
                } else {
                    *c = residue.to_lower_utf8_byte_amino()
                }
            });

        // match scores
        prf.emission_scores[Profile::MATCH_IDX]
            .iter_mut()
            .zip(mat_emissions)
            .skip(1)
            .for_each(|(scores, probs)| {
                scores
                    .iter_mut()
                    .zip(probs.iter())
                    .enumerate()
                    .take(Profile::MAX_ALPHABET_SIZE)
                    .for_each(|(i, (s, p))| {
                        // score is the (natural) log probability ratio:
                        //    ln(emission / background)
                        *s = (*p as f64 / AMINO_BACKGROUND_FREQUENCIES[i] as f64).ln() as f32
                    });

                // compute the ambiguity character
                // scores using the expected score
                let ambig = prf.alphabet.ambiguity_map();
                prf.alphabet
                    .ambiguous_iter()
                    .map(|a| (a as usize, &ambig[a]))
                    .filter(|(_, canon)| !canon.is_empty())
                    .for_each(|(a, canon)| {
                        let (res, denom) = canon.iter().map(|&c| c as usize).fold(
                            (0.0, 0.0),
                            |(res, denom), c| {
                                let weight = AMINO_BACKGROUND_FREQUENCIES[c];
                                (res + scores[c] * weight, denom + weight)
                            },
                        );
                        scores[a] = res / denom;
                    });
            });

        // insert scores
        prf.emission_scores[Profile::INSERT_IDX]
            .iter_mut()
            .skip(1)
            // insert at position M should be impossible,
            .take(prf.length - 1)
            // setting insert scores to 0 corresponds to insertion
            // emissions being equal to background probabilities
            //    ** because ln(P/P) = ln(1) = 0
            .for_each(|scores| scores.iter_mut().for_each(|s| *s = 0.0));

        Ok(prf)
    }
}

impl Profile {
    pub const MAX_ALPHABET_SIZE: usize = 20;
    pub const MAX_DEGENERATE_ALPHABET_SIZE: usize = 29;
    pub const GAP_INDEX: usize = 20;
    // non-residue character: "*"
    pub const NON_RESIDUE_IDX: usize = 27;
    // missing data character: "~"
    pub const MISSING_DATA_IDX: usize = 28;

    pub const MATCH_IDX: usize = 0;
    pub const INSERT_IDX: usize = 1;

    // special state indices
    pub const NUM_SPECIAL_STATES: usize = 5;
    pub const N_IDX: usize = 0;
    pub const B_IDX: usize = 1;
    pub const E_IDX: usize = 2;
    pub const C_IDX: usize = 3;
    pub const J_IDX: usize = 4;

    pub const SPECIAL_STATE_IDX_TO_NAME: [&'static str; 5] = ["N", "B", "E", "C", "J"];

    // special transition indices
    pub const SPECIAL_LOOP_IDX: usize = 0;
    pub const SPECIAL_MOVE_IDX: usize = 1;

    /// The number of allowed state transitions under the model.
    pub const NUM_STATE_TRANSITIONS: usize = 8;

    pub fn from_blosum_62_and_seq(seq: &Sequence) -> anyhow::Result<Self> {
        use blosum62::*;
        use Transition::*;

        let mut prf = ProfileBuilder::default();
        prf.length(seq.length);
        prf.name(seq.name.clone());

        // the transition probabilities are uniform
        let mut trans_buf = [0.0; 7];
        trans_buf[MM as usize] = 1.0 - 2.0 * P_OPEN;
        trans_buf[MI as usize] = P_OPEN;
        trans_buf[MD as usize] = P_OPEN;
        trans_buf[IM as usize] = 1.0 - P_EXTEND;
        trans_buf[II as usize] = P_EXTEND;
        trans_buf[DM as usize] = 1.0 - P_EXTEND;
        trans_buf[DD as usize] = P_EXTEND;

        // since we skip match emissions
        // at 0, throw these on first
        prf.transition(0, trans_buf);

        let mut emit_buf = [0.0f32; 20];

        seq.digital_bytes
            .iter()
            .enumerate()
            .skip(1)
            .take(seq.length)
            .for_each(|(pos, &residue)| {
                prf.transition(pos, trans_buf);
                emit_buf
                    .iter_mut()
                    .zip(CONDITIONAL_PROB[residue as usize])
                    .take(Alphabet::Amino.canonical_size())
                    .for_each(|(e, p)| *e = p);

                prf.mat_emission(pos, emit_buf);
            });

        // except for a few modifications to the last position
        trans_buf[MM as usize] = 1.0 - P_OPEN;
        trans_buf[MD as usize] = 0.0;
        trans_buf[DM as usize] = 1.0;
        trans_buf[DD as usize] = 0.0;
        prf.transition(seq.length, trans_buf);

        // lambda
        let total_entropy = prf
            .mat_emissions
            .iter()
            .skip(1)
            .map(|e| e.unwrap())
            .map(|probs| {
                probs
                    .iter()
                    .zip(AMINO_BACKGROUND_FREQUENCIES)
                    .map(|(p, f)| p * (p / f).log2())
                    .sum::<f32>()
            })
            .sum::<f32>();

        // note: we could get also get a pretty good approximation for
        //       lambda by simply using the product of the length of the
        //       sequence and the mean relative entropy of blosum62
        prf.fwd_lambda(std::f32::consts::LN_2 + 1.44 / total_entropy);

        // tau
        // note: this all comes from some experimentation fitting ~200k
        //       sequences to learn a closed form solution for tau
        const KL_BLOSUM62: f32 = 0.589;
        const MAGIC_TAU_1: f32 = 0.662;
        const MAGIC_TAU_2: f32 = 1.317;
        const MAGIC_TAU_3: f32 = 2.147;
        let tau = -((-MAGIC_TAU_1 * KL_BLOSUM62) + (MAGIC_TAU_2 * (seq.length as f32).ln())
            - MAGIC_TAU_3);

        prf.fwd_tau(tau);

        prf.build()
    }

    pub fn relative_entropy(&self) -> f32 {
        let probs: Vec<Vec<f32>> = self.emission_scores[Self::MATCH_IDX]
            .iter()
            // skip profile index 0
            .skip(1)
            // skip profile index L + 1
            .take(self.length)
            .map(|scores| {
                scores
                    .iter()
                    .take(20)
                    .enumerate()
                    .map(|(idx, score)| score.exp() * AMINO_BACKGROUND_FREQUENCIES[idx])
                    .collect::<Vec<f32>>()
            })
            .collect();

        mean_relative_entropy(&probs, &AMINO_BACKGROUND_FREQUENCIES)
    }

    pub fn adjust_mean_relative_entropy(&mut self, target_mre: f32) -> anyhow::Result<f32> {
        const LOWER_PROB_LIMIT: f32 = 1e-3;
        const TARGET_TOLERANCE: f32 = 1e-3;
        const WEIGHT_START: f32 = 0.1;
        const WEIGHT_STEP: f32 = 2.0;
        const MAX_ITER: usize = 100;

        let start_mre = self.relative_entropy();

        // double check that we aren't already at the target MRE
        if (start_mre - target_mre).abs() < TARGET_TOLERANCE {
            return Ok(start_mre);
        }

        let start_probs_by_pos: Vec<Vec<f32>> = self.emission_scores[Self::MATCH_IDX]
            .iter()
            .map(|scores| {
                scores
                    .iter()
                    .take(20)
                    .enumerate()
                    .map(|(idx, score)| score.exp() * AMINO_BACKGROUND_FREQUENCIES[idx])
                    .collect::<Vec<f32>>()
            })
            .collect();

        let mut new_probs_by_pos = start_probs_by_pos.clone();

        // for raising MRE, this will contain the maximum weights
        // that we can scale each position without setting any member
        // of the distribution below the LOWER_PROB_LIMIT
        let mut max_weights_by_pos: Vec<f32> = vec![0.0; self.length + 1];

        #[derive(Debug)]
        enum Mode {
            Raise,
            Lower,
        }

        let mode = if start_mre < target_mre {
            Mode::Raise
        } else {
            Mode::Lower
        };

        let (mut lower_bound, mut upper_bound) = match mode {
            Mode::Raise => {
                max_weights_by_pos
                    .iter_mut()
                    .enumerate()
                    // skip profile index 0
                    .skip(1)
                    .try_for_each(|(pos, weight)| {
                        *weight = start_probs_by_pos[pos]
                            .iter()
                            .enumerate()
                            .take(Profile::MAX_ALPHABET_SIZE)
                            .map(|(r, p)| {
                                // this comes from:
                                //   P'_a = (P_a + W * P_b) / (1 + W)
                                //   for P'_a = <LOWER_PROB_LIMIT>
                                (LOWER_PROB_LIMIT - p)
                                    / (AMINO_BACKGROUND_FREQUENCIES[r] - LOWER_PROB_LIMIT)
                            })
                            // **note: we're taking the max of (mostly) negative weights
                            .max_by(|a, b| a.partial_cmp(b).expect("NaN encountered"))
                            .ok_or(anyhow!("empty emission probabilities"))?;

                        anyhow::Ok(())
                    })
                    .context("failed to produce max weights")?;

                // the lower bound is the min of the max weights
                let lower_bound = max_weights_by_pos
                    .iter()
                    // skip profile index 0
                    .skip(1)
                    .min_by(|a, b| a.total_cmp(b))
                    .ok_or(anyhow!("empty weights"))?;

                Ok((*lower_bound, 0.0f32))
            }
            Mode::Lower => {
                let mut last_weight = 0.0;
                let mut weight = WEIGHT_START;
                let mut result = Err(anyhow!("exceeded max iterations during linear search"));
                for _ in 0..MAX_ITER {
                    new_probs_by_pos
                        .iter_mut()
                        .zip(&start_probs_by_pos)
                        // skip model position 0
                        .skip(1)
                        .take(self.length)
                        .for_each(|(new_probs, start_probs)| {
                            new_probs
                                .iter_mut()
                                .zip(start_probs)
                                .enumerate()
                                // only take core emission probs
                                .take(Profile::MAX_ALPHABET_SIZE)
                                // take a weighted mean
                                .for_each(|(residue_idx, (p_new, p_start))| {
                                    *p_new = (p_start
                                        + weight * AMINO_BACKGROUND_FREQUENCIES[residue_idx])
                                        / (1.0 + weight)
                                });
                        });

                    let current_mre = mean_relative_entropy(
                        &new_probs_by_pos[1..],
                        &AMINO_BACKGROUND_FREQUENCIES,
                    );

                    if current_mre < target_mre {
                        result = Ok((last_weight, weight));
                        break;
                    }

                    last_weight = weight;
                    weight *= WEIGHT_STEP;
                }
                result
            }
        }
        .with_context(|| {
            format!(
                "failed to produce weight bounds for MRE tuning: {}, {:4.3}->{:4.3}",
                match mode {
                    Mode::Raise => "raising",
                    Mode::Lower => "lowering",
                },
                start_mre,
                target_mre,
            )
        })?;

        let mut current_mre = start_mre;
        // subtle:
        //   by starting the weight at the lower bound instead of
        //   the mid point, we can easily catch the case of not
        //   being able to reach the target MRE when lowering
        let mut weight = lower_bound;

        for iter in 0..MAX_ITER {
            (1..=self.length).for_each(|pos| {
                let new_probs = &mut new_probs_by_pos[pos];
                let start_probs = &start_probs_by_pos[pos];
                let clamped_weight = weight.max(max_weights_by_pos[pos]);
                new_probs
                    .iter_mut()
                    .zip(start_probs)
                    .enumerate()
                    // only take core emission probs
                    .take(Profile::MAX_ALPHABET_SIZE)
                    .for_each(|(residue_idx, (p_new, p_start))| {
                        *p_new = (p_start
                            + clamped_weight * AMINO_BACKGROUND_FREQUENCIES[residue_idx])
                            / (1.0 + clamped_weight)
                    });
            });

            current_mre = mean_relative_entropy(
                &new_probs_by_pos[1..=self.length],
                &AMINO_BACKGROUND_FREQUENCIES,
            );

            let ordering = if (current_mre - target_mre).abs() < TARGET_TOLERANCE {
                Ordering::Equal
            } else {
                current_mre.total_cmp(&target_mre)
            };

            match ordering {
                Ordering::Less => upper_bound = weight,
                Ordering::Greater => lower_bound = weight,
                Ordering::Equal => break,
            }

            if lower_bound == upper_bound && iter == 0 {
                break;
            }

            weight = (lower_bound + upper_bound) / 2.0;
        }

        self.emission_scores[Self::MATCH_IDX]
            .iter_mut()
            .zip(new_probs_by_pos)
            // skip model position 0
            .skip(1)
            // skip model position L + 1
            .take(self.length)
            .for_each(|(scores, probs)| {
                scores
                    .iter_mut()
                    .zip(probs)
                    .enumerate()
                    .take(Profile::MAX_ALPHABET_SIZE)
                    .for_each(|(idx, (s, p))| *s = (p / AMINO_BACKGROUND_FREQUENCIES[idx]).ln())
            });

        Ok(current_mre)
    }

    pub fn calibrate_tau(&mut self, n: usize, target_length: usize, tail_probability: f32) {
        self.configure_for_target_length(target_length);

        let mut row_bounds = RowBounds::new(target_length);
        row_bounds.fill_rectangle(1, 1, target_length, self.length);
        let mut forward_matrix = DpMatrixSparse::new(target_length, self.length, &row_bounds);

        // first we are going to generate n random
        // sequences drawn from the background
        // distribution and compute their forward scores
        let mut scores = vec![0.0; n];
        let mut rng = Pcg64::seed_from_u64(0);
        (0..n).for_each(|seq_idx| {
            forward_matrix.reuse(target_length, self.length, &row_bounds);
            let seq = Sequence::random_amino(target_length, &mut rng);

            // **NOTE: HMMER uses multi-hit mode to calibrate
            //         Tau, but we are using uni-hit mode
            let forward_score_nats = forward(self, &seq, &mut forward_matrix, &row_bounds);

            let null_score_nats = null_one_score(target_length);

            scores[seq_idx] = (forward_score_nats - null_score_nats).to_bits().value();
        });

        /// This is some black magic from a textbook:
        ///   Statistical Models and Methods for
        ///   Lifetime Data by Joseph F. Lawless
        ///   
        /// From HMMER:
        ///   Equation 4.1.6 from [Lawless82], pg. 143, and
        ///   its first derivative with respect to lambda,
        ///   for finding the ML fit to Gumbel lambda parameter.
        ///   This equation gives a result of zero for the maximum
        ///   likelihood lambda.
        fn lawless416(samples: &[f32], lambda: f32) -> (f32, f32) {
            // e_sum is the sum of e^(-lambda x_i)
            let mut e_sum = 0.0f32;
            // x_sum is the sum of x_i
            let mut x_sum = 0.0f32;
            // xe_sum is the sum of x_i * e^(-lambda x_i)
            let mut xe_sum = 0.0f32;
            // xe_sum is the sum of x_i^2 * e^(-lambda x_i)
            let mut xxe_sum = 0.0f32;

            samples.iter().for_each(|x| {
                e_sum += (-lambda * x).exp();
                x_sum += x;
                xe_sum += x * (-lambda * x).exp();
                xxe_sum += x.powi(2) * (-lambda * x).exp();
            });

            let fx = (1.0 / lambda) - (x_sum / samples.len() as f32) + (xe_sum / e_sum);
            let dfx = (xe_sum / e_sum).powi(2) - (xxe_sum / e_sum) - (1.0 / (lambda.powi(2)));

            (fx, dfx)
        }

        // now make an initial guess at lambda
        //
        // from hmmer:
        //   (Evans/Hastings/Peacock, Statistical Distributions, 2000, p.86)
        let sum: f32 = scores.iter().sum();
        let squared_sum: f32 = scores.iter().map(|s| s.powi(2)).sum();
        let sample_variance: f32 =
            (squared_sum - sum * sum / scores.len() as f32) / (scores.len() as f32 - 1.0);
        let mut gumbel_lambda: f32 = std::f32::consts::PI / (6.0 / sample_variance).sqrt();

        // now we do this Newton/Raphson root finding
        // thing until we have converged on lambda
        let tolerance: f32 = 1e-5;
        let mut newton_raphson_success = false;
        for _ in 0..=100 {
            let (fx, dfx) = lawless416(&scores, gumbel_lambda);
            if fx.abs() < tolerance {
                newton_raphson_success = true;
                break;
            }
            gumbel_lambda -= fx / dfx;

            if gumbel_lambda <= 0.0 {
                gumbel_lambda = 0.001;
            }
        }

        if !newton_raphson_success {
            panic!("newton/raphson failed");
        }

        // this is apparently substituting into equation 4.1.5
        // from Lawless[82] to solve for the mu parameter
        let e_sum: f32 = scores.iter().map(|s| (-gumbel_lambda * s).exp()).sum();
        let gumbel_mu = -(e_sum / scores.len() as f32).ln() / gumbel_lambda;

        // now that we've fit the gumbel lambda, we are going to find our tau

        /// Calculates the inverse CDF for a Gumbel distribution
        /// with parameters <mu> and <lambda>. That is, returns
        /// the quantile <x> at which the CDF is <p>.
        fn gumbel_inverse_cdf(p: f32, mu: f32, lambda: f32) -> f32 {
            mu - ((-p.ln()).ln() / lambda)
        }

        // basically, there is (maybe?) no good method for fitting an exponential,
        // so instead we have fit a Gumbel that we are going to use to pick our tau
        //
        // from hmmer:
        //   Explanation of the eqn below: first find the x at which the Gumbel tail
        //   mass is predicted to be equal to tailp. Then back up from that x
        //   by log(tailp)/lambda to set the origin of the exponential tail to 1.0
        //   instead of tailp.
        self.fwd_tau = gumbel_inverse_cdf(1.0 - tail_probability, gumbel_mu, gumbel_lambda)
            + (tail_probability.ln() / self.fwd_lambda);
    }

    #[inline(always)]
    pub fn match_score(&self, alphabet_idx: usize, profile_idx: usize) -> f32 {
        self.emission_scores[Self::MATCH_IDX][profile_idx][alphabet_idx]
    }

    #[inline(always)]
    pub fn insert_score(&self, alphabet_idx: usize, profile_idx: usize) -> f32 {
        self.emission_scores[Self::INSERT_IDX][profile_idx][alphabet_idx]
    }

    #[inline(always)]
    pub fn transition_score(&self, transition_idx: usize, profile_idx: usize) -> f32 {
        self.core_transitions[profile_idx][transition_idx]
    }

    /// This is essentially a Kronecker delta function that returns 1.0 when the transition
    /// score is finite and f32::MIN when the transition score is -f32::INFINITY.
    ///
    /// Semantically, this means we are disallowing "impossible state paths" during posterior traceback.
    #[inline(always)]
    pub fn transition_score_delta(&self, transition_idx: usize, profile_idx: usize) -> f32 {
        if self.core_transitions[profile_idx][transition_idx] == -f32::INFINITY {
            f32::MIN_POSITIVE
        } else {
            1.0
        }
    }

    #[inline(always)]
    pub fn special_transition_score(&self, state_idx: usize, transition_idx: usize) -> f32 {
        self.special_transitions[state_idx][transition_idx]
    }

    /// This is essentially a Kronecker delta function that returns 1.0 when the transition
    /// score is finite and f32::MIN when the transition score is -f32::INFINITY.
    ///
    /// Semantically, this means we are disallowing "impossible state paths" during posterior traceback.
    #[inline(always)]
    pub fn special_transition_score_delta(&self, state_idx: usize, transition_idx: usize) -> f32 {
        if self.special_transitions[state_idx][transition_idx] == -f32::INFINITY {
            f32::MIN_POSITIVE
        } else {
            1.0
        }
    }

    pub fn generic_transition_score(
        &self,
        state_from: usize,
        idx_from: usize,
        state_to: usize,
        idx_to: usize,
    ) -> f32 {
        match state_from {
            Trace::S_STATE | Trace::T_STATE => 0.0,
            Trace::N_STATE => match state_to {
                Trace::B_STATE => {
                    self.special_transition_score(Profile::N_IDX, Profile::SPECIAL_MOVE_IDX)
                }
                Trace::N_STATE => {
                    self.special_transition_score(Profile::N_IDX, Profile::SPECIAL_LOOP_IDX)
                }
                _ => panic!(),
            },
            Trace::B_STATE => match state_to {
                Trace::M_STATE => self.transition_score(Transition::BM as usize, idx_to - 1),
                _ => panic!(),
            },
            Trace::M_STATE => match state_to {
                Trace::M_STATE => self.transition_score(Transition::MM as usize, idx_from),
                Trace::I_STATE => self.transition_score(Transition::MI as usize, idx_from),
                Trace::D_STATE => self.transition_score(Transition::MD as usize, idx_from),
                Trace::E_STATE => 0.0,
                _ => panic!(),
            },
            Trace::D_STATE => match state_to {
                Trace::M_STATE => self.transition_score(Transition::DM as usize, idx_from),
                Trace::D_STATE => self.transition_score(Transition::DD as usize, idx_from),
                Trace::E_STATE => 0.0,
                _ => panic!(),
            },
            Trace::I_STATE => match state_to {
                Trace::M_STATE => self.transition_score(Transition::IM as usize, idx_from),
                Trace::I_STATE => self.transition_score(Transition::II as usize, idx_from),
                _ => panic!(),
            },
            Trace::E_STATE => match state_to {
                Trace::C_STATE => {
                    self.special_transition_score(Profile::E_IDX, Profile::SPECIAL_MOVE_IDX)
                }
                Trace::J_STATE => {
                    self.special_transition_score(Profile::E_IDX, Profile::SPECIAL_LOOP_IDX)
                }
                _ => panic!(),
            },
            Trace::J_STATE => match state_to {
                Trace::B_STATE => {
                    self.special_transition_score(Profile::J_IDX, Profile::SPECIAL_MOVE_IDX)
                }
                Trace::J_STATE => {
                    self.special_transition_score(Profile::J_IDX, Profile::SPECIAL_LOOP_IDX)
                }
                _ => panic!(),
            },
            Trace::C_STATE => match state_to {
                Trace::T_STATE => {
                    self.special_transition_score(Profile::C_IDX, Profile::SPECIAL_MOVE_IDX)
                }
                Trace::C_STATE => {
                    self.special_transition_score(Profile::C_IDX, Profile::SPECIAL_LOOP_IDX)
                }
                _ => panic!(),
            },
            _ => panic!(),
        }
    }

    /// Sets the length of the current target sequence to which the profile will be aligned.
    ///
    /// This also adjusts the loop and move transition scores for the special states N, J, C.
    pub fn configure_for_target_length(&mut self, length: usize) {
        self.target_length = length;

        let move_probability: f32 =
            (2.0 + self.expected_j_uses) / (length as f32 + 2.0 + self.expected_j_uses);

        let loop_probability: f32 = 1.0 - move_probability;

        let loop_score = loop_probability.ln();
        let move_score = move_probability.ln();

        self.special_transitions[Profile::N_IDX][Profile::SPECIAL_LOOP_IDX] = loop_score;
        self.special_transitions[Profile::J_IDX][Profile::SPECIAL_LOOP_IDX] = loop_score;
        self.special_transitions[Profile::C_IDX][Profile::SPECIAL_LOOP_IDX] = loop_score;

        self.special_transitions[Profile::N_IDX][Profile::SPECIAL_MOVE_IDX] = move_score;
        self.special_transitions[Profile::J_IDX][Profile::SPECIAL_MOVE_IDX] = move_score;
        self.special_transitions[Profile::C_IDX][Profile::SPECIAL_MOVE_IDX] = move_score;

        // TODO: maybe temporary until I decide on whether
        //       or not to do a uniform entry distribution
        self.entry_transitions
            .iter_mut()
            .zip(self.core_transitions.iter())
            // skip position 0
            .skip(1)
            .for_each(|(a, b)| *a = b[Transition::BM as usize] + move_score);
    }
}

impl fmt::Debug for Profile {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        writeln!(f, "model length: {}", self.length)?;
        writeln!(f, "target length: {}", self.target_length)?;

        for i in 0..5 {
            writeln!(
                f,
                "{:8.4} {:8.4}",
                self.special_transitions[i][0], self.special_transitions[i][1]
            )?;
        }

        for i in 0..=self.length {
            writeln!(f, "{}", i)?;
            for residue in AMINO_ALPHABET_WITH_DEGENERATE {
                write!(f, "    {}   ", residue)?;
            }
            writeln!(f)?;

            for _ in 0..Profile::MAX_DEGENERATE_ALPHABET_SIZE {
                write!(f, "  ----- ")?;
            }
            writeln!(f)?;

            for j in 0..Profile::MAX_DEGENERATE_ALPHABET_SIZE {
                write!(f, "{:8.4} ", self.emission_scores[Self::MATCH_IDX][i][j])?;
            }
            writeln!(f)?;

            for j in 0..Profile::MAX_DEGENERATE_ALPHABET_SIZE {
                write!(f, "{:8.4} ", self.emission_scores[Self::INSERT_IDX][i][j])?;
            }
            writeln!(f)?;

            for t in 0..8 {
                write!(f, "{:8.4} ", self.core_transitions[i][t])?;
            }
            writeln!(f)?;
            writeln!(f)?;
        }

        Ok(())
    }
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

        let prf = Profile::from_blosum_62_and_seq(&seq)?;

        let composition_corrrect = [
            0.075741, 0.016550, 0.039336, 0.060816, 0.038081, 0.043406, 0.014495, 0.077001,
            0.061813, 0.144704, 0.026453, 0.033726, 0.040502, 0.031336, 0.048901, 0.063073,
            0.056054, 0.085715, 0.012350, 0.029948,
        ];

        todo!("fix this");
        // hmm.model
        //     .composition
        //     .iter()
        //     .zip(composition_corrrect)
        //     .for_each(|(a, b)| assert!((a - b).abs() < 1e-4));

        // let lambda_correct = 0.703323;
        // assert!((lambda_correct - hmm.stats.forward_lambda) < 1e-4,);

        Ok(())
    }

    #[test]
    fn test_calibrate_tau() -> anyhow::Result<()> {
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

        let mut prf = Profile::from_blosum_62_and_seq(&seq)?;

        prf.calibrate_tau(200, 100, 0.04);

        let correct_tau = -4.193134f32;
        let diff = (correct_tau - prf.fwd_tau).abs();

        // ***
        // run "cargo test -- --nocapture" to print the diff
        // ***
        println!("\x1b[31m{} | {}\x1b[0m", correct_tau, prf.fwd_tau);
        println!("\x1b[31mdifference: {}\x1b[0m", diff);
        assert!(diff <= 0.01);

        Ok(())
    }
}
