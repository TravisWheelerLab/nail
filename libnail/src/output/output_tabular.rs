use crate::align::structs::Alignment;

use anyhow::Result;
use std::io::Write;

const NUM_OUTPUT_COLUMNS: usize = 8;

const COLUMN_HEADERS: [&str; NUM_OUTPUT_COLUMNS] = [
    "target name",
    "profile name",
    "target start",
    "target end",
    "profile start",
    "profile end",
    "bit score",
    "e-value",
];

pub fn write_tabular_output(alignments: &Vec<Alignment>, out: &mut impl Write) -> Result<()> {
    let mut column_widths: [usize; NUM_OUTPUT_COLUMNS] = COLUMN_HEADERS.map(|s| s.len());
    for alignment in alignments {
        column_widths[0] = column_widths[0].max(alignment.target_name.len());
        column_widths[1] = column_widths[1].max(alignment.profile_name.len());
        column_widths[2] = column_widths[2].max(alignment.target_start.to_string().len());
        column_widths[3] = column_widths[3].max(alignment.target_end.to_string().len());
        column_widths[4] = column_widths[4].max(alignment.profile_start.to_string().len());
        column_widths[5] = column_widths[5].max(alignment.profile_end.to_string().len());
        column_widths[6] = column_widths[6].max(alignment.bit_score.to_string().len());
        column_widths[7] = column_widths[6].max(format!("{:1.1e}", alignment.evalue).len());
    }

    writeln!(
        out,
        "# {:w0$} {:w1$} {:w2$} {:w3$} {:w4$} {:w5$} {:w6$} {:w7$}",
        COLUMN_HEADERS[0],
        COLUMN_HEADERS[1],
        COLUMN_HEADERS[2],
        COLUMN_HEADERS[3],
        COLUMN_HEADERS[4],
        COLUMN_HEADERS[5],
        COLUMN_HEADERS[6],
        COLUMN_HEADERS[7],
        // subtract 2 from the first column width for the "# " prefix
        w0 = column_widths[0] - 2,
        w1 = column_widths[1],
        w2 = column_widths[2],
        w3 = column_widths[3],
        w4 = column_widths[4],
        w5 = column_widths[5],
        w6 = column_widths[6],
        w7 = column_widths[7],
    )?;

    writeln!(
        out,
        "#{} {} {} {} {} {} {} {}",
        // subtract 1 from the first width for the "#" prefix
        "-".repeat(column_widths[0] - 1),
        "-".repeat(column_widths[1]),
        "-".repeat(column_widths[2]),
        "-".repeat(column_widths[3]),
        "-".repeat(column_widths[4]),
        "-".repeat(column_widths[5]),
        "-".repeat(column_widths[6]),
        "-".repeat(column_widths[7]),
    )?;

    for alignment in alignments {
        writeln!(
            out,
            "{:w0$} {:w1$} {:w2$} {:w3$} {:w4$} {:w5$} {:w6$.2} {:w7$.1e}",
            alignment.target_name,
            alignment.profile_name,
            alignment.target_start,
            alignment.target_end,
            alignment.profile_start,
            alignment.profile_end,
            alignment.bit_score,
            alignment.evalue,
            w0 = column_widths[0],
            w1 = column_widths[1],
            w2 = column_widths[2],
            w3 = column_widths[3],
            w4 = column_widths[4],
            w5 = column_widths[5],
            w6 = column_widths[6],
            w7 = column_widths[7],
        )?;
    }

    Ok(())
}
