use crate::align::structs::Alignment;

use anyhow::Context;

trait ExtractString {
    fn extract_string(&self) -> String;
    fn empty_string(&self) -> String {
        "-".to_string()
    }
}

impl ExtractString for Option<String> {
    fn extract_string(&self) -> String {
        match self {
            Some(string) => string.to_string(),
            None => self.empty_string(),
        }
    }
}

impl ExtractString for Option<f32> {
    fn extract_string(&self) -> String {
        match self {
            Some(f) => format!("{f:.1}"),
            None => self.empty_string(),
        }
    }
}

impl ExtractString for Option<f64> {
    fn extract_string(&self) -> String {
        match self {
            Some(f) => format!("{f:.1e}"),
            None => self.empty_string(),
        }
    }
}

#[derive(Debug, Clone, Copy)]
pub enum Field {
    Target,
    Query,
    TargetStart,
    TargetEnd,
    QueryStart,
    QueryEnd,
    Score,
    CompBias,
    Evalue,
    CellFrac,
}

impl Field {
    pub fn extract_from(&self, alignment: &Alignment) -> String {
        match self {
            Field::Target => alignment.target_name.extract_string(),
            Field::Query => alignment.profile_name.extract_string(),
            Field::TargetStart => alignment.target_start.to_string(),
            Field::TargetEnd => alignment.target_end.to_string(),
            Field::QueryStart => alignment.profile_start.to_string(),
            Field::QueryEnd => alignment.profile_end.to_string(),
            Field::Score => alignment.score_bits.extract_string(),
            Field::CompBias => alignment.null_two.extract_string(),
            Field::Evalue => alignment.e_value.extract_string(),
            Field::CellFrac => alignment.cell_fraction.extract_string(),
        }
    }
}

#[derive(Clone)]
pub struct TableFormat {
    pub fields: Vec<Field>,
    pub labels: Vec<Vec<String>>,
    pub min_widths: Vec<usize>,
    pub widths: Vec<usize>,
}

impl TableFormat {
    pub fn new(fields: &[Field]) -> anyhow::Result<Self> {
        let mut labels = vec![];
        let mut min_widths = vec![];
        let mut widths = vec![];

        // this regex matches CamelCaseWords
        let label_regex =
            regex::Regex::new(r"[A-Z][a-z]*").context("failed to build field label regex")?;

        // this closure extracts the words & minimum column width for a field
        let label_fn = |field: &Field| -> anyhow::Result<(Vec<_>, usize), anyhow::Error> {
            // the Debug string for an enum produces the variant name
            let field_name = format!("{:?}", field);
            //
            // grab each word and its length in the variant name
            let (label_words, lengths): (Vec<_>, Vec<_>) = label_regex
                .find_iter(&field_name)
                .map(|m| (m.as_str().to_string().to_lowercase(), m.len()))
                .unzip();

            // the length of the longest word
            // is the min width of the column
            let min_width = *lengths
                .iter()
                .max()
                .context("failed to produce max field label width")?;
            Ok((label_words, min_width))
        };

        // we need to process the first field differently
        // because it needs to have at least +2 to
        // its minimum width to accomodate the "# " prefix
        let (mut label_words, mut min_width) = label_fn(&fields[0])?;
        labels.push(label_words);
        widths.push(min_width + 2);
        min_widths.push(min_width);

        for field in fields.iter().skip(1) {
            (label_words, min_width) = label_fn(field)?;
            labels.push(label_words);
            widths.push(min_width);
            min_widths.push(min_width);
        }

        Ok(Self {
            fields: fields.to_vec(),
            labels,
            min_widths,
            widths,
        })
    }

    pub fn update_widths(&mut self, alignment: &Alignment) {
        self.fields.iter().enumerate().for_each(|(idx, field)| {
            let width = field.extract_from(alignment).len();
            self.widths[idx] = self.widths[idx].max(width);
        });
    }

    pub fn reset_widths(&mut self) {
        self.widths
            .iter_mut()
            .zip(self.min_widths.iter())
            .for_each(|(width, min_width)| *width = *min_width);
    }

    pub fn header(&self) -> anyhow::Result<String> {
        // the number of rows in the header is
        // the max number of words in a field
        let num_rows = self
            // each entry in labels is
            // a vector of label words
            .labels
            .iter()
            .map(|l| l.len())
            .max()
            .context("field headers are empty")?;

        let mut header_row_strings: Vec<String> = vec!["# ".to_string(); num_rows + 1];

        // this function appends the field labels to the header
        let header_append_fn =
            |words: &Vec<String>, width: usize, header_row_strings: &mut Vec<String>| {
                let offset = num_rows - words.len();
                let mut words_padded = vec![""; offset];
                words.iter().for_each(|w| words_padded.push(w));

                words_padded.iter().enumerate().for_each(|(row, token)| {
                    let row_string = &mut header_row_strings[row];
                    *row_string = format!("{row_string}{:width$} ", token, width = width);
                });
                let last_row_string = header_row_strings.last_mut().unwrap();
                *last_row_string = format!("{last_row_string}{} ", "-".repeat(width));
            };

        // the first column gets -2 to it's width to account for the "# "
        header_append_fn(&self.labels[0], self.widths[0] - 2, &mut header_row_strings);

        self.labels
            .iter()
            // skip the first column
            .skip(1)
            .zip(self.widths.iter().skip(1))
            .for_each(|(words, &width)| {
                header_append_fn(words, width, &mut header_row_strings);
            });

        Ok(header_row_strings.join("\n"))
    }
}
