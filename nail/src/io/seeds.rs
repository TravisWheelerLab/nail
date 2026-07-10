use std::{
    collections::HashSet,
    fs::File,
    io::{BufRead, BufReader},
    path::Path,
};

use anyhow::bail;
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
                prf: tokens[0].to_string(),
                seq: tokens[1].to_string(),
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
        let mut prev = String::new();
        for line in buf.lines() {
            let line = line?;
            let tokens = line.split_whitespace().collect::<Vec<_>>();

            if tokens.len() != 8 {
                bail!("unexpected line format in mmseqs align tsv");
            }

            let prf = tokens[0].to_string();

            if prf != prev {
                if !seen.insert(prf.clone()) {
                    bail!("unexpected sorting of mmseqs align tsv");
                }
                prev = prf.clone();
                cnt = 0;
            }

            if cnt >= max {
                continue;
            }

            seeds.push(Seed {
                prf: prf,
                seq: tokens[1].to_string(),
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

#[cfg(test)]
mod tests {}
