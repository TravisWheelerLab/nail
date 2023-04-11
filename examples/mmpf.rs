use clap::Parser;
use nale::alphabet::AMINO_ALPHABET;
use std::fs::File;
use std::io::Read;
use std::path::PathBuf;

#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
pub struct Args {
    /// mmseqs profile database
    profile: String,
}

fn decode_floats_from_byte_chunk(byte_chunk: &[u8]) -> Option<Vec<f64>> {
    // every "row" in the mmseqs profile database
    // format should have 23 bytes
    if byte_chunk.len() != 23 {
        return None;
    }

    let mut floats: Vec<f64> = vec![];
    // the probabilities are stored in the first 20 bytes
    for byte in &byte_chunk[0..20] {
        // bias: 2^(k-1) - 1
        //       k = 3 (the number of bits in the exponent)
        // exponent with bias subtracted
        let exponent = (byte >> 5) as i32 - 3;
        // println!("{exponent}");
        // the mantissa has 5 bits
        let mantissa: u16 = (((byte << 3) >> 3) ^ 32) as u16;

        let decimal_point = 5 - exponent;
        // the start is basically "how many bits are we going to look at"
        // i.e. (the number of bits to the left) + (the number bits of bits to the right)
        let start = 6.max(decimal_point);

        let mut num = 0.0;

        // iterate across the bits from right to left
        // i.e. bit_idx = 0 is the rightmost bit
        for bit_idx in 0..start {
            // get the value of the bit
            let bit = (mantissa & (1 << bit_idx) != 0) as i32;
            let power = if bit_idx <= decimal_point {
                // if we are to the right of the decimal
                // point, we take a negative power
                -(decimal_point - bit_idx)
            } else {
                // if we are to the left of the decimal
                // point, we take a positive power
                bit_idx - decimal_point
            };

            // if the bit is set at the current bit_idx, we
            // want to include it as a component in the sum
            if bit == 1 {
                num += 2.0_f64.powi(power);
            }
        }
        floats.push(num);
    }
    Some(floats)
}

fn extract_query_character_from_byte_chunk(byte_chunk: &[u8]) -> Option<String> {
    if byte_chunk.len() == 23 {
        Some(String::from(AMINO_ALPHABET[byte_chunk[20] as usize]))
    } else {
        None
    }
}

fn extract_consensus_character_from_byte_chunk(byte_chunk: &[u8]) -> Option<String> {
    if byte_chunk.len() == 23 {
        Some(String::from(AMINO_ALPHABET[byte_chunk[21] as usize]))
    } else {
        None
    }
}

fn main() {
    let args = Args::parse();
    let path = PathBuf::from(&args.profile);
    let mut file = File::open(path).unwrap();
    let mut buffer = Vec::new();
    file.read_to_end(&mut buffer).unwrap();

    let mut query_seq = String::from("");
    let mut consensus_seq = String::from("");
    for byte_chunk in buffer.chunks(23) {
        // if let Some(floats) = decode_floats_from_byte_chunk(byte_chunk) {
        //     for val in floats {
        //         print!("{:8.3} ", val)
        //     }
        //     println!();
        // }
        // if let Some(s) = extract_query_character_from_byte_chunk(byte_chunk) {
        //     query_seq.push_str(&s);
        // }
        if let Some(s) = extract_consensus_character_from_byte_chunk(byte_chunk) {
            consensus_seq.push_str(&s);
        }
    }

    println!("{consensus_seq}");
    println!("{query_seq}");
}
