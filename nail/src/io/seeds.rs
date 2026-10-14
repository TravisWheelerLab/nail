use std::{
    collections::HashSet,
    fs::File,
    io::{BufRead, BufReader},
    path::Path,
};

use anyhow::{bail, Context};
use libnail::align::structs::Seed;

pub struct Seeds {
    pub seeds: Vec<Seed>,
}

impl Seeds {
    pub fn from_path<P: AsRef<Path>>(path: P, max_seeds: Option<usize>) -> anyhow::Result<Self> {
        let file = BufReader::new(File::open(path)?);

        if let Some(max) = max_seeds {
            Self::new_max(file, max)
        } else {
            Self::new(file)
        }
    }

    fn new<R: BufRead>(buf: R) -> anyhow::Result<Self> {
        let mut seeds = vec![];
        for line in buf.lines() {
            let line = line?;
            let tokens = line.split_whitespace().collect::<Vec<_>>();
            seeds.push(Seed {
                prf_pos: position(tokens[0])?,
                seq_pos: position(tokens[1])?,
                seq_start: tokens[4].parse()?,
                seq_end: tokens[5].parse()?,
                prf_start: tokens[2].parse()?,
                prf_end: tokens[3].parse()?,
                score: tokens[6].parse()?,
                e_value: tokens[7].parse()?,
            })
        }
        Ok(Self { seeds })
    }

    fn new_max<R: BufRead>(buf: R, max: usize) -> anyhow::Result<Self> {
        let mut seeds = vec![];
        let mut seen = HashSet::new();

        let mut cnt = 0;
        let mut prev = None;
        for line in buf.lines() {
            let line = line?;
            let tokens = line.split_whitespace().collect::<Vec<_>>();

            if tokens.len() != 8 {
                bail!("unexpected line format in mmseqs align tsv");
            }

            let prf = position(tokens[0])?;

            if prev != Some(prf) {
                if !seen.insert(prf) {
                    bail!("unexpected sorting of mmseqs align tsv");
                }
                prev = Some(prf);
                cnt = 0;
            }

            if cnt >= max {
                continue;
            }

            seeds.push(Seed {
                prf_pos: prf,
                seq_pos: position(tokens[1])?,
                seq_start: tokens[4].parse()?,
                seq_end: tokens[5].parse()?,
                prf_start: tokens[2].parse()?,
                prf_end: tokens[3].parse()?,
                score: tokens[6].parse()?,
                e_value: tokens[7].parse()?,
            });

            cnt += 1;
        }

        Ok(Self { seeds })
    }
}

fn position(token: &str) -> anyhow::Result<usize> {
    token
        .parse()
        .with_context(|| format!("seed field \"{token}\" is not a record position"))
}

#[cfg(test)]
mod tests {}
