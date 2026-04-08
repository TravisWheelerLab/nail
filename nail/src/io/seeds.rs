use std::{
    fs::File,
    io::{BufRead, BufReader},
    path::Path,
};

use libnail::align::structs::Seed;

pub struct Seeds {
    pub seeds: Vec<Seed>,
}

impl Seeds {
    pub fn from_path<P: AsRef<Path>>(path: P) -> anyhow::Result<Self> {
        let file = BufReader::new(File::open(path)?);

        let mut seeds = vec![];
        for line in file.lines() {
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
}

#[cfg(test)]
mod tests {}
