use phf::phf_map;

pub const UTF8_SPACE: u8 = 32;
pub const UTF8_STAR: u8 = 42;
pub const UTF8_PLUS: u8 = 43;
pub const UTF8_DASH: u8 = 45;
pub const UTF8_DOT: u8 = 46;
pub const UTF8_PIPE: u8 = 124;

/// maps from \<usize\> -> \<UTF8 value for usize\>
pub const UTF8_NUMERIC: [u8; 11] = [48, 49, 50, 51, 52, 53, 54, 55, 56, 57, UTF8_STAR];

pub const AMINO_ALPHABET: [&str; 20] = [
    "A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W",
    "Y",
];

pub const AMINO_ALPHABET_WITH_DEGENERATE: [&str; 29] = [
    "A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W",
    "Y", "-", "B", "J", "Z", "O", "U", "X", "*", "~",
];

pub const UTF8_TO_DIGITAL_AMINO: phf::Map<u8, u8> = phf_map! {
    // upper case
    65u8 => 0,    // A
    67u8 => 1,    // C
    68u8 => 2,    // D
    69u8 => 3,    // E
    70u8 => 4,    // F
    71u8 => 5,    // G
    72u8 => 6,    // H
    73u8 => 7,    // I
    75u8 => 8,    // K
    76u8 => 9,    // L
    77u8 => 10,   // M
    78u8 => 11,   // N
    80u8 => 12,   // P
    81u8 => 13,   // Q
    82u8 => 14,   // R
    83u8 => 15,   // S
    84u8 => 16,   // T
    86u8 => 17,   // V
    87u8 => 18,   // W
    89u8 => 19,   // Y
    // lower case
    97u8 => 0,    // a
    99u8 => 1,    // c
    100u8 => 2,   // d
    101u8 => 3,   // e
    102u8 => 4,   // f
    103u8 => 5,   // g
    104u8 => 6,   // h
    105u8 => 7,   // i
    107u8 => 8,   // k
    108u8 => 9,   // l
    109u8 => 10,  // m
    110u8 => 11,  // n
    112u8 => 12,  // p
    113u8 => 13,  // q
    114u8 => 14,  // r
    115u8 => 15,  // s
    116u8 => 16,  // t
    118u8 => 17,  // v
    119u8 => 18,  // w
    121u8 => 19,  // y
    // degenerate characters
    79u8 => 20,   // O
    85u8 => 21,   // U
    88u8 => 22,   // X
    66u8 => 23,   // B
    90u8 => 24,   // Z
    74u8 => 25,   // J
    111u8 => 20,  // o
    117u8 => 21,  // u
    120u8 => 22,  // x
    98u8 => 23,   // b
    122u8 => 24,  // z
    106u8 => 25,  // j
};

pub const AMINO_INVERSE_MAP: phf::Map<u8, u8> = phf_map! {
    0u8  => 65,   // A
    1u8  => 67,   // C
    2u8  => 68,   // D
    3u8  => 69,   // E
    4u8  => 70,   // F
    5u8  => 71,   // G
    6u8  => 72,   // H
    7u8  => 73,   // I
    8u8  => 75,   // K
    9u8  => 76,   // L
    10u8 => 77,   // M
    11u8 => 78,   // N
    12u8 => 80,   // P
    13u8 => 81,   // Q
    14u8 => 82,   // R
    15u8 => 83,   // S
    16u8 => 84,   // T
    17u8 => 86,   // V
    18u8 => 87,   // W
    19u8 => 89,   // Y
    // end base alphabet
    20u8 => 79,   // O
    21u8 => 85,   // U
    22u8 => 88,   // X
    23u8 => 66,   // B
    24u8 => 90,   // Z
    25u8 => 74,   // J
    45u8 => 45,   // -
    46u8 => 46,   // .
    32u8 => 32,   // space
    255u8 => 32,  // space
};

pub const AMINO_INVERSE_MAP_LOWER: phf::Map<u8, u8> = phf_map! {
    0u8  => 97,   // a
    1u8  => 99,   // c
    2u8  => 100,   // d
    3u8  => 101,   // e
    4u8  => 102,   // f
    5u8  => 103,   // g
    6u8  => 104,   // h
    7u8  => 105,   // i
    8u8  => 107,   // k
    9u8  => 108,   // l
    10u8 => 109,   // m
    11u8 => 110,   // n
    12u8 => 112,   // p
    13u8 => 113,   // q
    14u8 => 114,   // r
    15u8 => 115,   // s
    16u8 => 116,   // t
    17u8 => 118,   // v
    18u8 => 119,   // w
    19u8 => 121,   // y
    // end base alphabet
    111u8 => 79,   // o
    117u8 => 85,   // u
    120u8 => 88,   // x
    98u8  => 66,   // b
    122u8 => 90,   // z
    106u8 => 74,   // j
    45u8 => 45,    // -
    46u8 => 46,    // .
    32u8 => 32,    // space
    255u8 => 32,   // space
};

// TODO: where did these come from?
pub const AMINO_BACKGROUND_FREQUENCIES: [f32; 20] = [
    0.0787945, // A
    0.0151600, // C
    0.0535222, // D
    0.0668298, // E
    0.0397062, // F
    0.0695071, // G
    0.0229198, // H
    0.0590092, // I
    0.0594422, // K
    0.0963728, // L
    0.0237718, // M
    0.0414386, // N
    0.0482904, // P
    0.0395639, // Q
    0.0540978, // R
    0.0683364, // S
    0.0540687, // T
    0.0673417, // V
    0.0114135, // W
    0.0304133, // Y
];
