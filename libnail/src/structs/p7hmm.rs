use anyhow::{anyhow, bail, Context};

use super::profile::{ProfileBuilder, Transition};
use super::Profile;

/// The transitions in the order a p7hmm file lists them.
const P7HMM_TRANSITIONS: [Transition; 7] = [
    Transition::MM,
    Transition::MI,
    Transition::MD,
    Transition::IM,
    Transition::II,
    Transition::DM,
    Transition::DD,
];

impl Profile {
    /// Build a profile from the bytes of one p7hmm record.
    pub fn from_p7hmm(buf: &[u8]) -> anyhow::Result<Profile> {
        enum ParseState {
            Header,
            Model,
        }

        #[derive(Debug)]
        enum ModelState {
            Comp,
            Mat,
            Ins,
            Trans,
        }

        let mut emission_buf = [0.0; Profile::MAX_ALPHABET_SIZE];
        let mut trans_buf = [0.0; Profile::NUM_STATE_TRANSITIONS - 1];
        let mut swap_buf = [0.0; Profile::NUM_STATE_TRANSITIONS - 1];
        let mut state = ParseState::Header;
        let mut model_state = ModelState::Comp;

        let mut prf = ProfileBuilder::default();
        let mut pos = 0;
        for line in buf.split(|b| *b == b'\n').filter(|b| !b.is_empty()) {
            match state {
                ParseState::Header => {
                    let line = std::str::from_utf8(line)?.trim();

                    // a header line is a keyword, then its value
                    let (keyword, value) = match line.split_once(char::is_whitespace) {
                        Some((keyword, value)) => (keyword, value.trim()),
                        None => (line, ""),
                    };

                    match keyword {
                        "NAME" => {
                            prf.name(value.to_string());
                        }
                        "ACC" => {
                            prf.accession(value.to_string());
                        }
                        "LENG" => {
                            prf.length(value.parse().with_context(|| {
                                format!("failed to parse HMM length: \"{}\"", value)
                            })?);
                        }
                        "ALPH" => {
                            if value != "amino" {
                                bail!("non-amino HMM alphabet found: \"{}\"", value);
                            };
                        }
                        "STATS" => {
                            let stats: Vec<&str> = value.split_whitespace().collect();

                            if stats[1] == "FORWARD" {
                                prf.fwd_tau(stats[2].parse().with_context(|| {
                                    format!("failed to parse HMM Forward tau: \"{}\"", stats[2])
                                })?);
                                prf.fwd_lambda(stats[3].parse().with_context(|| {
                                    format!("failed to parse HMM Forward lambda: \"{}\"", stats[3])
                                })?);
                            }
                        }
                        "HMM" => state = ParseState::Model,
                        _ => {
                            // ignore most flags for now
                        }
                    }
                }
                ParseState::Model => {
                    match model_state {
                        ModelState::Comp => {
                            if first_word(line)? == "COMPO" {
                                // note: we don't actually use the
                                //       composition in the p7HMM
                                //  why: HMMER doesn't
                                model_state = ModelState::Ins;
                            }
                        }
                        ModelState::Mat => {
                            let pos_pos = match first_non_whitespace_pos(line) {
                                Some(pos) => pos,
                                // line was empty
                                None => continue,
                            };
                            let pos_str = match word_from(line, pos_pos)? {
                                "//" => {
                                    break;
                                }
                                s => s,
                            };

                            let new_pos = pos_str.parse::<usize>()?;

                            if new_pos != pos + 1 {
                                bail!("non-linear HMM positions: {pos}->{new_pos}");
                            } else {
                                pos = new_pos
                            }

                            parse_p7hmm_floats_into::<{ Profile::MAX_ALPHABET_SIZE }>(
                                &line[pos_pos + pos_str.len()..],
                                &mut emission_buf,
                            )?;

                            prf.mat_emission(pos, emission_buf);
                            model_state = ModelState::Ins;
                        }
                        ModelState::Ins => {
                            // note: we don't actually use the
                            //       insert emissions in the p7HMM
                            //  why: HMMER doesn't
                            model_state = ModelState::Trans;
                        }
                        ModelState::Trans => {
                            parse_p7hmm_floats_into::<{ Profile::NUM_STATE_TRANSITIONS - 1 }>(
                                line,
                                &mut trans_buf,
                            )?;

                            // map p7HMM transition order to Profile order
                            for (p7hmm_idx, transition) in P7HMM_TRANSITIONS.into_iter().enumerate()
                            {
                                swap_buf[transition as usize] = trans_buf[p7hmm_idx];
                            }

                            prf.transition(pos, swap_buf);

                            model_state = ModelState::Mat;
                        }
                    }
                }
            }
        }

        prf.build()
    }
}

fn negexp(x: f32) -> f32 {
    (-x).exp()
}

fn parse_p7hmm_floats_into<const N: usize>(
    floats: &[u8],
    out: &mut [f32; N],
) -> anyhow::Result<()> {
    let mut tokens = floats
        .split(|b| b.is_ascii_whitespace())
        .filter(|token| !token.is_empty());

    for item in out.iter_mut() {
        *item = match tokens.next() {
            Some(b"*") => 0.0,
            Some(token) => negexp(parse_p7hmm_float(token)?),
            None => bail!("expected {N} floats, got fewer"),
        };
    }

    Ok(())
}

/// Parse one float of a p7hmm file, giving the same `f32` as `str::parse`.
fn parse_p7hmm_float(token: &[u8]) -> anyhow::Result<f32> {
    // what: read the digits as an integer n and count the
    //       digits after the point, d; the value is n / 10^d
    // why: with at most 7 digits, n < 2^24 and 10^d are both
    //      exact in f32, and IEEE division rounds correctly,
    //      so the quotient is the f32 nearest the decimal,
    //      which is what str::parse returns
    const MAX_DIGITS: usize = 7;
    const POWERS_OF_TEN: [f32; MAX_DIGITS + 1] = [1e0, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7];

    let mut n: u32 = 0;
    let mut digits = 0;
    let mut point_seen = false;
    let mut fraction_digits = 0;

    for &byte in token {
        if byte.is_ascii_digit() && digits < MAX_DIGITS {
            n = n * 10 + (byte - b'0') as u32;
            digits += 1;
            if point_seen {
                fraction_digits += 1;
            }
        } else if byte == b'.' && !point_seen && digits > 0 {
            point_seen = true;
        } else {
            // a sign, an exponent, a leading or second point,
            // or an eighth digit
            return parse_with_str(token);
        }
    }

    if digits == 0 {
        return parse_with_str(token);
    }

    Ok(n as f32 / POWERS_OF_TEN[fraction_digits])
}

fn parse_with_str(token: &[u8]) -> anyhow::Result<f32> {
    let token = std::str::from_utf8(token)?;
    token
        .parse()
        .with_context(|| format!("failed to parse HMM float: \"{token}\""))
}

fn word_from(buf: &[u8], pos: usize) -> anyhow::Result<&str> {
    let end = buf[pos..]
        .iter()
        .position(|b| b.is_ascii_whitespace())
        .map(|n| pos + n)
        .unwrap_or(buf.len());

    Ok(std::str::from_utf8(&buf[pos..end])?)
}

fn first_non_whitespace_pos(buf: &[u8]) -> Option<usize> {
    buf.iter().position(|b| !b.is_ascii_whitespace())
}

fn first_word(buf: &[u8]) -> anyhow::Result<&str> {
    word_from(
        buf,
        first_non_whitespace_pos(buf).ok_or(anyhow!("no words found in buffer"))?,
    )
}

#[cfg(test)]
mod tests {
    use std::{
        fs::File,
        io::{Cursor, Read},
        path::Path,
    };

    use super::*;

    fn assert_parses_like_str_parse(token: &str) {
        match (parse_p7hmm_float(token.as_bytes()), token.parse::<f32>()) {
            (Ok(ours), Ok(std)) => assert_eq!(ours.to_bits(), std.to_bits(), "{token}"),
            (Err(_), Err(_)) => {}
            (ours, std) => panic!("{token}: ours {ours:?}, str::parse {std:?}"),
        }
    }

    #[test]
    fn parse_p7hmm_float_matches_str_parse_on_every_five_decimal_value() {
        // p7hmm writes its floats with five decimals, so this
        // covers every value from 0.00000 to 99.99999
        for n in 0..10_000_000 {
            assert_parses_like_str_parse(&format!("{}.{:05}", n / 100_000, n % 100_000));
        }
    }

    #[test]
    fn parse_p7hmm_float_matches_str_parse_at_each_point_position() {
        for decimals in 0..=8 {
            let scale = 10u32.pow(decimals);
            for n in (0..100_000_000).step_by(997) {
                let token = match decimals {
                    0 => format!("{n}"),
                    _ => format!(
                        "{}.{:0width$}",
                        n / scale,
                        n % scale,
                        width = decimals as usize
                    ),
                };
                assert_parses_like_str_parse(&token);
            }
        }
    }

    #[test]
    fn parse_p7hmm_float_matches_str_parse_on_other_tokens() {
        for token in [
            "0",
            "1.",
            ".5",
            "0.00000",
            "9999999",
            "99999999",
            "16777217",
            "1.2345678",
            "+1.0",
            "-1.0",
            "1e-3",
            "inf",
            "",
            ".",
            "-",
            "1.2.3",
            "abc",
        ] {
            assert_parses_like_str_parse(token);
        }
    }

    #[test]
    fn test_profile_serialize() -> anyhow::Result<()> {
        let mut bytes = vec![];
        let mut file =
            File::open(Path::new(env!("CARGO_MANIFEST_DIR")).join("../fixtures/query.hmm"))?;
        file.read_to_end(&mut bytes)?;
        let prf = Profile::from_p7hmm(&bytes)?;

        let mut buf = vec![];

        prf.serialize(&mut buf)?;

        let de = Profile::deserialize(Cursor::new(buf))?;

        assert_eq!(prf.name, de.name);
        assert_eq!(prf.accession, de.accession);
        assert_eq!(prf.consensus_seq_bytes_utf8, de.consensus_seq_bytes_utf8);
        assert_eq!(prf.core_transitions, de.core_transitions);
        assert_eq!(prf.emission_scores[0], de.emission_scores[0]);
        Ok(())
    }
}
