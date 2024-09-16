use std::io::{stdout, Write};

use libnail::{
    align::structs::Alignment,
    output::output_tabular::{Field, TableFormat},
};

use crate::{args::OutputArgs, util::PathBufExt};

pub const DEFAULT_COLUMNS: [Field; 10] = [
    Field::Target,
    Field::Query,
    Field::TargetStart,
    Field::TargetEnd,
    Field::QueryStart,
    Field::QueryEnd,
    Field::Score,
    Field::CompBias,
    Field::Evalue,
    Field::CellFrac,
];

pub enum HeaderStatus {
    Unwritten,
    Written,
}

pub struct OutputStep {
    alignment_writer: Box<dyn Write + Send>,
    table_writer: Box<dyn Write + Send>,
    table_format: TableFormat,
    header_status: HeaderStatus,
}

impl OutputStep {
    pub fn new(args: &OutputArgs) -> anyhow::Result<Self> {
        Ok(Self {
            alignment_writer: match &args.ali_results_path {
                Some(path) => Box::new(path.open(true)?),
                None => Box::new(stdout()),
            },
            table_writer: Box::new(args.tbl_results_path.open(true)?),
            table_format: TableFormat::new(&DEFAULT_COLUMNS)?,
            header_status: HeaderStatus::Unwritten,
        })
    }

    pub fn write(&mut self, alignments: &mut [Alignment]) -> anyhow::Result<()> {
        self.table_format.reset_widths();
        self.table_format.update_widths(alignments);

        alignments.sort_by(|a, b| a.scores.e_value.partial_cmp(&b.scores.e_value).unwrap());

        if let HeaderStatus::Unwritten = self.header_status {
            let header = TableFormat::header(&self.table_format)?;
            writeln!(self.table_writer, "{header}")?;
            self.header_status = HeaderStatus::Written;
        }

        alignments.iter().for_each(|ali| {
            writeln!(
                self.table_writer,
                "{}",
                ali.tab_string_formatted(&self.table_format)
            )
            .expect("failed to write tabular output");

            writeln!(self.alignment_writer, "{}", ali.ali_string())
                .expect("failed to write alignment output");
        });

        Ok(())
    }
}
