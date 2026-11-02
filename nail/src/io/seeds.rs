use std::{
    collections::HashSet,
    fs::File,
    io::{BufRead, BufReader, Write},
    path::Path,
    str::FromStr,
};

use anyhow::{bail, Context};
use libnail::align::structs::Seed;

pub struct Seeds {
    pub seeds: Vec<Seed>,
}

/// Turns a seeds file's query or target field into a record position.
pub type FieldToPosition<'a> = &'a dyn Fn(&str) -> anyhow::Result<usize>;

impl Seeds {
    pub fn from_path<P: AsRef<Path>>(
        path: P,
        max_seeds: Option<usize>,
        query_position: FieldToPosition,
        target_position: FieldToPosition,
    ) -> anyhow::Result<Self> {
        let path = path.as_ref();
        let file = File::open(path)
            .with_context(|| format!("failed to open seeds file: {}", path.display()))?;
        let file = BufReader::new(file);

        let seeds = match max_seeds {
            Some(max) => Self::new_max(file, max, query_position, target_position),
            None => Self::new(file, query_position, target_position),
        };
        seeds.with_context(|| format!("failed to read seeds file: {}", path.display()))
    }

    fn new<R: BufRead>(
        buf: R,
        query_position: FieldToPosition,
        target_position: FieldToPosition,
    ) -> anyhow::Result<Self> {
        let mut seeds = vec![];
        for (line_idx, line) in buf.lines().enumerate() {
            let line = line?;
            if line.trim().is_empty() {
                continue;
            }
            seeds.push(parse_seed(
                line_idx + 1,
                &line,
                query_position,
                target_position,
            )?);
        }
        Ok(Self { seeds })
    }

    fn new_max<R: BufRead>(
        buf: R,
        max: usize,
        query_position: FieldToPosition,
        target_position: FieldToPosition,
    ) -> anyhow::Result<Self> {
        let mut seeds = vec![];
        let mut seen = HashSet::new();

        let mut cnt = 0;
        let mut prev = None;
        for (line_idx, line) in buf.lines().enumerate() {
            let line = line?;
            if line.trim().is_empty() {
                continue;
            }
            let seed = parse_seed(line_idx + 1, &line, query_position, target_position)?;

            if prev != Some(seed.prf_pos) {
                if !seen.insert(seed.prf_pos) {
                    bail!(
                        "line {}: seeds for query {} are not contiguous",
                        line_idx + 1,
                        seed.prf_pos
                    );
                }
                prev = Some(seed.prf_pos);
                cnt = 0;
            }

            if cnt >= max {
                continue;
            }

            seeds.push(seed);
            cnt += 1;
        }

        Ok(Self { seeds })
    }
}

fn parse_seed(
    line_number: usize,
    line: &str,
    query_position: FieldToPosition,
    target_position: FieldToPosition,
) -> anyhow::Result<Seed> {
    let fields = line.split_whitespace().collect::<Vec<_>>();
    if fields.len() != 8 {
        bail!("line {line_number} has {} fields, expected 8", fields.len());
    }

    Ok(Seed {
        prf_pos: query_position(fields[0])?,
        seq_pos: target_position(fields[1])?,
        prf_start: parse_number(line_number, &fields, 2)?,
        prf_end: parse_number(line_number, &fields, 3)?,
        seq_start: parse_number(line_number, &fields, 4)?,
        seq_end: parse_number(line_number, &fields, 5)?,
        score: parse_number(line_number, &fields, 6)?,
        e_value: parse_number(line_number, &fields, 7)?,
    })
}

fn parse_number<T>(line_number: usize, fields: &[&str], field_idx: usize) -> anyhow::Result<T>
where
    T: FromStr,
    T::Err: std::error::Error + Send + Sync + 'static,
{
    let field = fields[field_idx];
    field.parse().with_context(|| {
        format!(
            "line {line_number}, field {}: \"{field}\" is not a number",
            field_idx + 1
        )
    })
}

/// Copy a seeds file that holds positions, naming each query and target.
//
// the fields after the first two are copied as they are,
// so the copy differs from its source only in those two
pub fn write_named_seeds<R: BufRead, W: Write>(
    positional: R,
    mut named: W,
    query_name: &dyn Fn(usize) -> anyhow::Result<String>,
    target_name: &dyn Fn(usize) -> anyhow::Result<String>,
) -> anyhow::Result<()> {
    for line in positional.lines() {
        let line = line?;
        let mut fields = line.split('\t');

        let query = fields.next().context("a seeds line is empty")?;
        let target = fields.next().context("a seeds line has no target field")?;
        write!(
            named,
            "{}\t{}",
            query_name(parse_position(query)?)?,
            target_name(parse_position(target)?)?
        )?;

        for field in fields {
            write!(named, "\t{field}")?;
        }
        writeln!(named)?;
    }

    named.flush()?;
    Ok(())
}

/// Read a seeds file field that holds a record position.
pub fn parse_position(token: &str) -> anyhow::Result<usize> {
    token
        .parse()
        .with_context(|| format!("seed field \"{token}\" is not a record position"))
}

#[cfg(test)]
mod tests {
    use super::*;

    // two seeds in the format mmseqs convertalis writes:
    // query, target, then six fields nail copies through
    const POSITIONAL: &str = "0\t2\t1\t259\t54\t302\t189\t3.660E-55\n\
                              1\t0\t3\t120\t10\t130\t95\t1.281E-20\n";

    const QUERY_NAMES: [&str; 2] = ["7tm_1", "7tm_2"];
    const TARGET_NAMES: [&str; 3] = ["O08858", "F1MV99", "A0T2N3"];

    fn query_name(n: usize) -> anyhow::Result<String> {
        Ok(QUERY_NAMES[n].to_string())
    }

    fn target_name(n: usize) -> anyhow::Result<String> {
        Ok(TARGET_NAMES[n].to_string())
    }

    fn query_position(name: &str) -> anyhow::Result<usize> {
        QUERY_NAMES
            .iter()
            .position(|&q| q == name)
            .with_context(|| format!("no query {name}"))
    }

    fn target_position(name: &str) -> anyhow::Result<usize> {
        TARGET_NAMES
            .iter()
            .position(|&t| t == name)
            .with_context(|| format!("no target {name}"))
    }

    #[test]
    fn named_seeds_change_only_the_first_two_fields() -> anyhow::Result<()> {
        let mut named = vec![];
        write_named_seeds(POSITIONAL.as_bytes(), &mut named, &query_name, &target_name)?;

        assert_eq!(
            String::from_utf8(named)?,
            "7tm_1\tA0T2N3\t1\t259\t54\t302\t189\t3.660E-55\n\
             7tm_2\tO08858\t3\t120\t10\t130\t95\t1.281E-20\n"
        );

        Ok(())
    }

    #[test]
    fn named_seeds_read_back_as_the_positional_ones() -> anyhow::Result<()> {
        let mut named = vec![];
        write_named_seeds(POSITIONAL.as_bytes(), &mut named, &query_name, &target_name)?;

        let from_positions = Seeds::new(POSITIONAL.as_bytes(), &parse_position, &parse_position)?;
        let from_names = Seeds::new(named.as_slice(), &query_position, &target_position)?;

        assert_eq!(from_names.seeds, from_positions.seeds);
        assert_eq!(from_positions.seeds[0].prf_pos, 0);
        assert_eq!(from_positions.seeds[0].seq_pos, 2);

        Ok(())
    }

    #[test]
    fn a_short_line_names_the_line_and_its_field_count() {
        let seeds = "0\t2\t1\t259\t54\t302\t189\t3.660E-55\n\
                     0\t1\t3\t120\n";

        let Err(err) = Seeds::new(seeds.as_bytes(), &parse_position, &parse_position) else {
            panic!("a four-field line was accepted");
        };
        assert_eq!(err.to_string(), "line 2 has 4 fields, expected 8");
    }

    #[test]
    fn a_non_number_names_the_line_and_field() {
        let seeds = "0\t2\t1\t259\tfifty\t302\t189\t3.660E-55\n";

        let Err(err) = Seeds::new(seeds.as_bytes(), &parse_position, &parse_position) else {
            panic!("a word was read as a number");
        };
        assert_eq!(
            err.to_string(),
            "line 1, field 5: \"fifty\" is not a number"
        );
    }

    #[test]
    fn blank_lines_are_skipped() -> anyhow::Result<()> {
        let seeds = "\n\
                     0\t2\t1\t259\t54\t302\t189\t3.660E-55\n\
                     \n\
                     1\t0\t3\t120\t10\t130\t95\t1.281E-20\n\
                     \n";

        let seeds = Seeds::new(seeds.as_bytes(), &parse_position, &parse_position)?;
        assert_eq!(seeds.seeds.len(), 2);
        assert_eq!(seeds.seeds[1].prf_pos, 1);
        Ok(())
    }

    #[test]
    fn a_name_in_a_positional_seeds_file_fails() {
        let named = "7tm_1\tA0T2N3\t1\t259\t54\t302\t189\t3.660E-55\n";

        let Err(err) = Seeds::new(named.as_bytes(), &parse_position, &parse_position) else {
            panic!("a name was read as a position");
        };
        assert!(err
            .to_string()
            .contains("\"7tm_1\" is not a record position"));
    }
}
