use crate::align::structs::Alignment;

use anyhow::Result;
use std::io::Write;

pub fn write_standard_output(alignments: &Vec<Alignment>, out: &mut impl Write) -> Result<()> {
    for alignment in alignments {
        writeln!(out, "{}", &alignment.ali_string())?
    }
    Ok(())
}
