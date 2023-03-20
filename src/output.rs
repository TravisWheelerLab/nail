use crate::structs::Alignment;
use anyhow::Result;
use std::io::Write;

pub fn write_tabular_output(alignments: &Vec<Alignment>, out: &mut impl Write) -> Result<()> {
    for alignment in alignments {
        writeln!(
            out,
            "{} {} {} {}",
            alignment.target_name,
            alignment.profile_name,
            alignment.target_start,
            alignment.target_end
        )?;
    }

    Ok(())
}
