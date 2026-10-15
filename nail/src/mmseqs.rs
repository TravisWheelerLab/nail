use self::consts::*;

use std::{
    fs::File,
    io::{BufRead, BufReader, Write},
    ops::Range,
    path::{Path, PathBuf},
    process::Command,
};

use libnail::{align::Nats, alphabet::UTF8_TO_DIGITAL_AMINO, structs::Profile};
use libsail::source::{FileSource, Source};

use crate::{
    args::SearchArgs,
    io::Fasta,
    util::{CommandExt, PathExt},
};

use anyhow::Context;
use anyhow::{anyhow, bail};
use rayon::prelude::*;
use regex::Regex;

pub mod consts {
    pub const AMINO_DBTYPE: &[u8] = &[0, 0, 0, 0];
    pub const PROFILE_DBTYPE: &[u8] = &[2, 0, 0, 0];
    pub const ALIGN_DBTYPE: &[u8] = &[5, 0, 0, 0];
    pub const PREFILTER_DBTYPE: &[u8] = &[7, 0, 0, 0];
    pub const GENERIC_DBTYPE: &[u8] = &[12, 0, 0, 0];
}

#[derive(Debug, PartialEq)]
pub enum ByteBuffer {
    Complete(Vec<u8>),
    Partial(Vec<u8>, usize),
    Empty,
}

#[derive(Clone, Copy)]
struct Descriptor {
    file_idx: usize,
    offset: u64,
    length: u64,
}

pub struct PrefilterDb {
    /// The path the db was opened from.
    path: PathBuf,
    /// The split files, in suffix order.
    sources: Vec<FileSource>,
    descriptors: Vec<Descriptor>,
    checkpoints: Vec<u64>,
}

impl PrefilterDb {
    pub fn from_path<P>(path: P) -> anyhow::Result<Self>
    where
        P: AsRef<Path>,
    {
        let path = path.as_ref();
        let dir = path.parent().context("path has no parent dir")?;
        let name = path
            .file_stem()
            .context("no file stem")?
            .to_str()
            .context("invalid utf8")?;

        // ---

        let index_path = path.with_extension("index");
        let index_reader =
            BufReader::new(File::open(&index_path).with_context(|| format!("{:?}", index_path))?);

        let mut offsets = vec![];
        let mut lengths = vec![];
        for line in index_reader.lines() {
            let line = line?;

            let tokens = line
                .split('\t')
                .map(|s| s.parse::<u64>().expect("failed to parse index line"))
                .collect::<Vec<_>>();

            offsets.push(tokens[1]);
            // -1 because the index lengths include a null byte
            lengths.push(tokens[2] - 1);
        }

        offsets.shrink_to_fit();
        lengths.shrink_to_fit();

        // ---

        let re = Regex::new(&format!(r"^{}\.\d+$", regex::escape(name)))?;
        let mut paths = std::fs::read_dir(dir)?
            .filter_map(Result::ok)
            .map(|e| e.path())
            .filter(|p| {
                p.file_name()
                    .and_then(|s| s.to_str())
                    .map(|s| re.is_match(s))
                    .unwrap_or(false)
            })
            .collect::<Vec<_>>();

        paths.sort_by_key(|p| {
            p.extension()
                .and_then(|e| e.to_str())
                .and_then(|s| s.parse::<usize>().ok())
                .expect("bad DB split suffix")
        });

        // if we have no paths here, that
        // means the prefilter is one file
        if paths.is_empty() {
            paths.push(path.to_path_buf())
        }

        let file_sizes: Vec<u64> = paths
            .iter()
            .map(|p| p.metadata().map(|m| m.len()))
            .collect::<Result<_, _>>()?;

        let prefix_sum = file_sizes
            .into_iter()
            .scan(0u64, |acc, x| {
                let cur = *acc;
                *acc += x;
                Some(cur)
            })
            .collect::<Vec<_>>();

        // --

        let descriptors = offsets
            .into_iter()
            .zip(lengths)
            .map(|(offset, length)| {
                // the last file starting at or before the offset:
                // an empty split file starts where the next one
                // does, and binary_search may return either for
                // a record at the start of the next file
                let file_idx = prefix_sum.partition_point(|&start| start <= offset) - 1;

                let relative_offset = offset - prefix_sum[file_idx];
                Descriptor {
                    file_idx,
                    offset: relative_offset,
                    length,
                }
            })
            .collect::<Vec<_>>();

        let sources = paths
            .iter()
            .map(|path| {
                FileSource::open(path).with_context(|| format!("failed to open {}", path.display()))
            })
            .collect::<anyhow::Result<Vec<_>>>()?;

        let checkpoints = vec![0; descriptors.len()];

        Ok(Self {
            path: path.to_path_buf(),
            sources,
            descriptors,
            checkpoints,
        })
    }

    pub fn get(&self, prf_idx: usize) -> anyhow::Result<Vec<u8>> {
        let desc = self.descriptors[prf_idx];
        let length = desc.length as usize;

        // the null byte is read with the record
        let mut bytes = self.sources[desc.file_idx]
            .range(desc.offset, length + 1)?
            .into_owned();
        if bytes[length] != 0 {
            bail!(
                "a record in {} does not end in a null byte",
                self.path.display()
            );
        }

        bytes.truncate(length);
        Ok(bytes)
    }

    pub fn next_n(&mut self, prf_idx: usize, n: usize) -> anyhow::Result<ByteBuffer> {
        let desc = self.descriptors[prf_idx];
        let checkpoint = &mut self.checkpoints[prf_idx];

        if n == 0 {
            return Ok(ByteBuffer::Empty);
        }

        let start = desc.offset + *checkpoint;
        let remaining = (desc.length - *checkpoint) as usize;

        // the null byte is checked when a take reaches the
        // end of the record, and not for every record in
        // from_path: that cost a disk read per record, and
        // 6.9 s on a cold 32 GB db through a memory map
        if remaining == 0 {
            check_null_byte(&self.sources, desc, &self.path)?;
            return Ok(ByteBuffer::Empty);
        }

        let source = &self.sources[desc.file_idx];
        let mut bytes = vec![];
        let mut chunk_len = SCAN_CHUNK;
        let mut record_cnt = 0;

        while bytes.len() < remaining {
            let scanned = bytes.len();
            let len = chunk_len.min(remaining - scanned);
            bytes.extend_from_slice(&source.range(start + scanned as u64, len)?);

            for pos in scanned..bytes.len() {
                if bytes[pos] == b'\n' {
                    record_cnt += 1;
                }

                if record_cnt == n {
                    *checkpoint += (pos + 1) as u64;
                    bytes.truncate(pos + 1);
                    return Ok(ByteBuffer::Complete(bytes));
                }
            }

            chunk_len *= 2;
        }

        check_null_byte(&self.sources, desc, &self.path)?;
        *checkpoint = desc.length;
        Ok(ByteBuffer::Partial(bytes, record_cnt))
    }
}

/// The length of the first range `next_n` reads; each later range doubles.
//
// a take reads less than twice the bytes it needs,
// plus one first range, so a small take reads little
const SCAN_CHUNK: usize = 1 << 12;

/// Fail unless the byte after a record is its null byte.
fn check_null_byte(sources: &[FileSource], desc: Descriptor, path: &Path) -> anyhow::Result<()> {
    let byte = sources[desc.file_idx].range(desc.offset + desc.length, 1)?;
    if byte[0] != 0 {
        bail!("a record in {} does not end in a null byte", path.display());
    }

    Ok(())
}

pub struct MmseqsDbPaths {
    pub query_db: PathBuf,
    pub target_db: PathBuf,
    pub prefilter_db: PathBuf,
    pub align_db: PathBuf,
    pub prog_dir: PathBuf,
}

impl MmseqsDbPaths {
    pub fn new(dir: impl AsRef<Path>) -> Self {
        let dir = dir.as_ref();

        Self {
            query_db: dir.join("query-db/qdb"),
            target_db: dir.join("target-db/tdb"),
            prefilter_db: dir.join("prefilter-db/pdb"),
            align_db: dir.join("align-db/adb"),
            prog_dir: dir.join("prog-seed/"),
        }
    }

    pub fn destroy(&self) -> anyhow::Result<()> {
        Self::remove_db(&self.query_db).context("failed to remove mmseqs query DB")?;
        Self::remove_db(&self.target_db).context("failed to remove mmseqs target DB")?;
        Self::remove_db(&self.prefilter_db).context("failed to remove mmseqs prefilter DB")?;
        Self::remove_db(&self.align_db).context("failed to remove mmseqs align DB")?;

        if self
            .prog_dir
            .try_exists()
            .with_context(|| format!("failed to check existence of: {:?}", self.prog_dir))?
        {
            let re = Regex::new(r"^\d+$")?;
            std::fs::read_dir(&self.prog_dir)
                .with_context(|| format!("failed to open dir: {:?}", &self.prog_dir))?
                .try_for_each(|entry| -> anyhow::Result<()> {
                    let entry = entry.with_context(|| {
                        format!("failed to access a file in dir: {:?}", &self.prog_dir)
                    })?;

                    let path = entry.path();

                    let ft = entry
                        .file_type()
                        .with_context(|| format!("failed to stat {path:?}"))?;

                    if ft.is_dir()
                        && path
                            .file_name()
                            .and_then(|s| s.to_str())
                            .is_some_and(|s| re.is_match(s))
                    {
                        std::fs::remove_dir_all(&path)
                            .with_context(|| format!("failed to remove {path:?}"))?;
                    }

                    Ok(())
                })?;
        }

        Ok(())
    }

    pub fn check(&self) -> anyhow::Result<()> {
        if Self::db_exists(&self.query_db)? {
            bail!("mmseqs query DB already exists");
        } else if Self::db_exists(&self.target_db)? {
            bail!("mmseqs target DB already exists");
        } else if Self::db_exists(&self.prefilter_db)? {
            bail!("mmseqs prefilter DB already exists");
        } else if Self::db_exists(&self.align_db)? {
            bail!("mmseqs align DB already exists");
        }
        Ok(())
    }

    fn db_exists<P: AsRef<Path>>(path: P) -> anyhow::Result<bool> {
        let path = path.as_ref();
        let db_dir = path.parent().context("DB path has no parent dir")?;

        if !db_dir
            .try_exists()
            .with_context(|| format!("failed to check existence of: {db_dir:?}"))?
        {
            return Ok(false);
        }

        let db_name = path
            .file_stem()
            .and_then(|s| s.to_str())
            .context("failed to get db name from path")?;

        let re = Regex::new(&format!(r"^{}.*", regex::escape(db_name)))?;

        std::fs::read_dir(db_dir)
            .with_context(|| format!("failed to open dir: {db_dir:?}"))?
            .try_fold(false, |_, entry| {
                let p = entry
                    .with_context(|| format!("failed to access a file in dir: {db_dir:?}"))?
                    .path();

                Ok(p.file_name()
                    .and_then(|s| s.to_str())
                    .is_some_and(|s| re.is_match(s)))
            })
    }

    fn remove_db<P: AsRef<Path>>(path: P) -> anyhow::Result<()> {
        let path = path.as_ref();
        let db_dir = path.parent().context("DB path has no parent dir")?;

        if !db_dir
            .try_exists()
            .with_context(|| format!("failed to check existence of: {db_dir:?}"))?
        {
            return Ok(());
        }

        if !db_dir.is_dir() {
            bail!("DB parent is not a directory");
        }

        if !db_dir
            .file_name()
            .and_then(|s| s.to_str())
            .is_some_and(|s| s.ends_with("-db"))
        {
            bail!("DB parent dir name is supposed to end with -db");
        }

        let db_name = path
            .file_stem()
            .and_then(|s| s.to_str())
            .context("failed to get db name from path")?;

        if db_name.is_empty() {
            bail!("empty DB prefix")
        }

        let re = Regex::new(&format!(
            r"^{}(\.\d+|\.index|\.dbtype|_h|_h\.index|_h\.dbtype)?$",
            regex::escape(db_name)
        ))?;

        std::fs::read_dir(db_dir)
            .with_context(|| format!("failed to open dir: {db_dir:?}"))?
            .try_for_each(|entry| -> anyhow::Result<()> {
                let entry =
                    entry.with_context(|| format!("failed to access a file in dir: {db_dir:?}"))?;

                let path = entry.path();

                let ft = entry
                    .file_type()
                    .with_context(|| format!("failed to stat {path:?}"))?;

                if ft.is_dir() {
                    return Ok(());
                }

                if path
                    .file_name()
                    .and_then(|s| s.to_str())
                    .is_some_and(|s| re.is_match(s))
                {
                    std::fs::remove_file(&path)
                        .with_context(|| format!("failed to remove {path:?}"))?;
                }

                Ok(())
            })?;

        Ok(())
    }
}

/// The number of records in one batch of a sequence db.
const SEQUENCE_DB_BATCH_LEN: usize = 16_384;

/// The bytes a run of records adds to a sequence db and its header db.
struct SequenceDbBatch {
    db: Vec<u8>,
    db_lengths: Vec<usize>,
    headers: Vec<u8>,
    header_lengths: Vec<usize>,
}

impl SequenceDbBatch {
    fn build(sequences: &Fasta, positions: Range<usize>) -> anyhow::Result<Self> {
        let mut batch = Self {
            db: vec![],
            db_lengths: Vec::with_capacity(positions.len()),
            headers: vec![],
            header_lengths: Vec::with_capacity(positions.len()),
        };

        for position in positions {
            let seq = sequences.get(position)?;
            batch.db.extend_from_slice(&seq.utf8_bytes[1..]);
            batch.db.extend_from_slice(&[b'\n', 0]);
            batch.db_lengths.push(seq.length + 2);

            // the header is the record's position, so the seeds
            // mmseqs reports name records by position
            let header = position.to_string();
            batch.headers.extend_from_slice(header.as_bytes());
            batch.headers.extend_from_slice(&[b'\n', 0]);
            batch.header_lengths.push(header.len() + 2);
        }

        Ok(batch)
    }
}

pub fn write_mmseqs_sequence_database(
    sequences: &Fasta,
    path: impl AsRef<Path>,
) -> anyhow::Result<()> {
    let db_path = path.as_ref().to_owned();
    let db_name = db_path.file_name().unwrap().to_str().unwrap();
    let db_index_path = db_path.with_file_name(format!("{db_name}.index"));
    let db_dbtype_path = db_path.with_file_name(format!("{db_name}.dbtype"));

    let header_path = db_path.with_file_name(format!("{db_name}_h"));
    let header_index_path = db_path.with_file_name(format!("{db_name}_h.index"));
    let header_dbtype_path = db_path.with_file_name(format!("{db_name}_h.dbtype"));

    let db_dir = db_path
        .parent()
        .context("failed to produce mmseqs sequence DB directory path")?;

    db_dir.create_dir()?;

    db_dbtype_path.open(true)?.write_all(AMINO_DBTYPE)?;
    header_dbtype_path.open(true)?.write_all(GENERIC_DBTYPE)?;

    let mut db = db_path.open(true)?;
    let mut db_index = db_index_path.open(true)?;
    let mut db_header = header_path.open(true)?;
    let mut db_header_index = header_index_path.open(true)?;

    write_sequence_db(
        sequences,
        SEQUENCE_DB_BATCH_LEN,
        &mut db,
        &mut db_index,
        &mut db_header,
        &mut db_header_index,
    )
}

fn write_sequence_db(
    sequences: &Fasta,
    batch_len: usize,
    db: &mut impl Write,
    db_index: &mut impl Write,
    headers: &mut impl Write,
    header_index: &mut impl Write,
) -> anyhow::Result<()> {
    let mut db_offset = 0usize;
    let mut header_offset = 0usize;
    let mut position = 0usize;

    let batch_starts: Vec<usize> = (0..sequences.len()).step_by(batch_len).collect();

    // a group of batches is built at a time, so memory
    // holds one group's bytes and not the whole db
    for group in batch_starts.chunks(rayon::current_num_threads()) {
        let batches = group
            .par_iter()
            .map(|&start| {
                let end = (start + batch_len).min(sequences.len());
                SequenceDbBatch::build(sequences, start..end)
            })
            .collect::<anyhow::Result<Vec<_>>>()?;

        for batch in batches {
            db.write_all(&batch.db)?;
            headers.write_all(&batch.headers)?;

            for (db_length, header_length) in batch.db_lengths.iter().zip(&batch.header_lengths) {
                writeln!(db_index, "{position}\t{db_offset}\t{db_length}")?;
                writeln!(header_index, "{position}\t{header_offset}\t{header_length}")?;

                db_offset += db_length;
                header_offset += header_length;
                position += 1;
            }
        }
    }

    Ok(())
}

pub fn write_mmseqs_profile_database(
    profiles: &[Profile],
    path: impl AsRef<Path>,
) -> anyhow::Result<()> {
    let db_path = path.as_ref().to_owned();
    let db_name = db_path.file_name().unwrap().to_str().unwrap();
    let db_index_path = db_path.with_file_name(format!("{db_name}.index"));
    let db_dbtype_path = db_path.with_file_name(format!("{db_name}.dbtype"));

    let header_path = db_path.with_file_name(format!("{db_name}_h"));
    let header_index_path = db_path.with_file_name(format!("{db_name}_h.index"));
    let header_dbtype_path = db_path.with_file_name(format!("{db_name}_h.dbtype"));

    let db_dir = db_path
        .parent()
        .context("failed to produce mmseqs profile DB directory path")?;

    db_dir.create_dir()?;

    db_dbtype_path.open(true)?.write_all(PROFILE_DBTYPE)?;
    header_dbtype_path.open(true)?.write_all(GENERIC_DBTYPE)?;

    let mut db = db_path.open(true)?;
    let mut db_index = db_index_path.open(true)?;
    let mut db_header = header_path.open(true)?;
    let mut db_header_index = header_index_path.open(true)?;

    let mut db_offset = 0usize;
    let mut header_offset = 0usize;

    for (prf_cnt, prf) in profiles.iter().enumerate() {
        for prf_idx in 1..=prf.length {
            for byte in (0..20)
                .map(|residue| Nats(prf.match_score(residue, prf_idx)))
                .map(|nats| nats.to_bits())
                .map(|bits| bits.value())
                // multiply the bits by 8.0 just because
                // that's what they do in mmseqs
                .map(|s| (s * 8.0) as i8)
                .map(|s| s as u8)
            {
                db.write_all(&[byte])?;
            }

            let consensus_byte_digital = *UTF8_TO_DIGITAL_AMINO
                .get(&prf.consensus_seq_bytes_utf8[prf_idx])
                .unwrap();

            // query sequence byte?
            // still not sure what this is, so we'll just put the consensus
            db.write_all(&[consensus_byte_digital])?;

            // consensus sequence byte
            db.write_all(&[consensus_byte_digital])?;

            // neff value
            db.write_all(&[0u8])?;

            // something
            db.write_all(&[0u8])?;
            // something
            db.write_all(&[0u8])?;
        }

        let db_byte_length = prf.length * 25;

        writeln!(db_index, "{}\t{}\t{}", prf_cnt, db_offset, db_byte_length,)?;

        // for some reason, the header has newlines and 0-byte separators?
        let header = prf_cnt.to_string();
        writeln!(db_header, "{header}")?;
        db_header.write_all(&[0u8])?;

        // +1 for the 0 byte, +1 for the newline
        let header_byte_length = header.len() + 2;

        writeln!(
            db_header_index,
            "{}\t{}\t{}",
            prf_cnt, header_offset, header_byte_length
        )?;

        db_offset += db_byte_length;
        header_offset += header_byte_length;
    }

    Ok(())
}

pub fn run_mmseqs_prefilter(
    query_db_path: impl AsRef<Path>,
    target_db_path: impl AsRef<Path>,
    prefilter_db_path: impl AsRef<Path>,
    // TODO: generics make this annoying
    score_mx_path: Option<PathBuf>,
    args: &SearchArgs,
) -> anyhow::Result<()> {
    let mut prefilter = Command::new(&args.mmseqs_path);

    let qdb = query_db_path.as_ref();
    let tdb = target_db_path.as_ref();
    let pdb = prefilter_db_path.as_ref();

    let pdb_dir = pdb
        .parent()
        .context("failed to produce mmseqs prefilter DB directory path")?;

    pdb_dir.create_dir()?;

    prefilter
        .arg("prefilter")
        .arg(qdb)
        .arg(tdb)
        .arg(pdb)
        .args(["--threads", &args.num_threads.to_string()])
        .args(["-k", &args.seed_args.k.to_string()])
        .args(["-s", &args.seed_args.s.to_string()])
        .args(["--max-seqs", &args.seed_args.max_seqs.to_string()]);

    if let Some(v) = args.seed_args.comp_bias_corr {
        prefilter.args(["--comp-bias-corr", &v.to_string()]);
    }

    if let Some(ref path) = score_mx_path {
        prefilter.arg("--sub-mat");
        prefilter.arg(path);
    };

    prefilter.run()?;

    Ok(())
}

pub fn run_mmseqs_align(
    query_db_path: impl AsRef<Path>,
    target_db_path: impl AsRef<Path>,
    prefilter_db_path: impl AsRef<Path>,
    align_db_path: impl AsRef<Path>,
    // TODO: generics make this annoying
    score_mx_path: Option<PathBuf>,
    args: &SearchArgs,
) -> anyhow::Result<()> {
    let effective_e_value = args.pipeline_args.seed_pvalue_threshold
        * args
            .expert_args
            .target_database_size
            .ok_or(anyhow!("no target database size"))? as f64;

    let qdb = query_db_path.as_ref();
    let tdb = target_db_path.as_ref();
    let pdb = prefilter_db_path.as_ref();
    let adb = align_db_path.as_ref().to_path_buf();

    let adb_dir = adb
        .parent()
        .context("failed to produce mmseqs align DB directory path")?;

    adb_dir.create_dir()?;

    let mut align = Command::new(&args.mmseqs_path);
    align
        .arg("align")
        .arg(qdb)
        .arg(tdb)
        .arg(pdb)
        .arg(&adb)
        .args(["--threads", &args.num_threads.to_string()])
        .args(["-e", &effective_e_value.to_string()])
        // the '-a' argument enables alignment backtraces in mmseqs2
        // it is required to get start positions for alignments
        .args(["-a", "1"]);

    if let Some(v) = args.seed_args.comp_bias_corr {
        align.args(["--comp-bias-corr", &v.to_string()]);
    }

    if let Some(ref path) = score_mx_path {
        align.arg("--sub-mat");
        align.arg(path);
    };

    align.run()?;

    Ok(())
}

pub fn run_mmseqs_convertalis(
    query_db_path: impl AsRef<Path>,
    target_db_path: impl AsRef<Path>,
    align_db_path: impl AsRef<Path>,
    align_tsv_path: impl AsRef<Path>,
    args: &SearchArgs,
) -> anyhow::Result<()> {
    let qdb = query_db_path.as_ref();
    let tdb = target_db_path.as_ref();
    let adb = align_db_path.as_ref().to_path_buf();
    let align_tsv = align_tsv_path.as_ref();

    Command::new(&args.mmseqs_path)
        .arg("convertalis")
        .arg(qdb)
        .arg(tdb)
        .arg(adb)
        .arg(align_tsv)
        .args(["--threads", &args.num_threads.to_string()])
        .args([
            "--format-output",
            "qheader,theader,qstart,qend,tstart,tend,bits,evalue",
        ])
        .run()?;

    Ok(())
}

#[cfg(test)]
mod tests {
    use std::path::Path;

    use super::*;

    // prefilter-db holds three records over three split
    // files. pdb.0 holds record 0, "a\nb\nc\n". pdb.1 is
    // empty, so it starts at offset 7, as pdb.2 does. pdb.2
    // holds record 1, "d\n", and record 2, which is empty
    fn split_prefilter_db() -> anyhow::Result<PrefilterDb> {
        let path = Path::new(env!("CARGO_MANIFEST_DIR")).join("../fixtures/prefilter-db/pdb");
        PrefilterDb::from_path(path)
    }

    #[test]
    fn prefilter_db_gets_records_across_split_files() -> anyhow::Result<()> {
        let db = split_prefilter_db()?;

        assert_eq!(db.get(0)?, b"a\nb\nc\n");
        assert_eq!(db.get(1)?, b"d\n");
        assert_eq!(db.get(2)?, b"");

        Ok(())
    }

    #[test]
    fn prefilter_db_next_n_takes_lines_from_a_checkpoint() -> anyhow::Result<()> {
        let mut db = split_prefilter_db()?;

        // a take of zero lines does not move the checkpoint
        assert_eq!(db.next_n(0, 0)?, ByteBuffer::Empty);
        assert_eq!(db.next_n(0, 2)?, ByteBuffer::Complete(b"a\nb\n".to_vec()));

        // one line is left, so the take is partial and counts it
        assert_eq!(db.next_n(0, 2)?, ByteBuffer::Partial(b"c\n".to_vec(), 1));
        assert_eq!(db.next_n(0, 2)?, ByteBuffer::Empty);

        // each record keeps its own checkpoint
        assert_eq!(db.next_n(1, 1)?, ByteBuffer::Complete(b"d\n".to_vec()));
        assert_eq!(db.next_n(1, 1)?, ByteBuffer::Empty);
        assert_eq!(db.next_n(2, 1)?, ByteBuffer::Empty);

        Ok(())
    }

    // the file is "a\n\0", and the index gives the one
    // record a length of 2 with its null byte. so the
    // record is "a", and the byte after it is "\n"
    fn unterminated_prefilter_db() -> anyhow::Result<PrefilterDb> {
        let path =
            Path::new(env!("CARGO_MANIFEST_DIR")).join("../fixtures/prefilter-db-unterminated/pdb");
        PrefilterDb::from_path(path)
    }

    #[test]
    fn prefilter_db_get_refuses_a_record_without_its_null_byte() -> anyhow::Result<()> {
        let db = unterminated_prefilter_db()?;

        let Err(err) = db.get(0) else {
            panic!("a record without its null byte was accepted");
        };
        assert!(err.to_string().contains("does not end in a null byte"));

        Ok(())
    }

    #[test]
    fn prefilter_db_next_n_refuses_a_record_without_its_null_byte() -> anyhow::Result<()> {
        let mut db = unterminated_prefilter_db()?;

        // "a" has no newline, so a take of one line reaches
        // the end of the record and checks the byte after it
        let Err(err) = db.next_n(0, 1) else {
            panic!("a record without its null byte was accepted");
        };
        assert!(err.to_string().contains("does not end in a null byte"));

        Ok(())
    }

    #[test]
    fn sequence_db_is_the_same_at_every_batch_length() -> anyhow::Result<()> {
        let path = Path::new(env!("CARGO_MANIFEST_DIR")).join("../fixtures/target.fa");
        let fasta = Fasta::from_path(&path)?;

        // the db as the writer built it before batching:
        // one record at a time, in file order. text, so a
        // failure prints lines and not bytes
        let mut expected_db = String::new();
        let mut expected_db_index = String::new();
        let mut expected_headers = String::new();
        let mut expected_header_index = String::new();
        let mut db_offset = 0;
        let mut header_offset = 0;
        for n in 0..fasta.len() {
            let seq = fasta.get(n)?;
            expected_db += std::str::from_utf8(&seq.utf8_bytes[1..])?;
            expected_db += "\n\0";
            expected_db_index += &format!("{n}\t{db_offset}\t{}\n", seq.length + 2);
            db_offset += seq.length + 2;

            let header = n.to_string();
            expected_headers += &format!("{header}\n\0");
            expected_header_index += &format!("{n}\t{header_offset}\t{}\n", header.len() + 2);
            header_offset += header.len() + 2;
        }

        // 1 is a batch per record, 7 leaves a last batch of
        // one record, 400 is one batch, 401 overruns the file
        for batch_len in [1, 7, 400, 401] {
            // at 3 threads a group is 3 batches, and the last
            // group is short at batch lengths 1 and 7
            for threads in [1, 3] {
                let pool = rayon::ThreadPoolBuilder::new()
                    .num_threads(threads)
                    .build()?;

                let mut db = vec![];
                let mut db_index = vec![];
                let mut headers = vec![];
                let mut header_index = vec![];
                pool.install(|| {
                    write_sequence_db(
                        &fasta,
                        batch_len,
                        &mut db,
                        &mut db_index,
                        &mut headers,
                        &mut header_index,
                    )
                })?;

                let at = format!("batch length {batch_len}, {threads} threads");
                assert_eq!(String::from_utf8(db)?, expected_db, "db at {at}");
                assert_eq!(
                    String::from_utf8(db_index)?,
                    expected_db_index,
                    ".index at {at}"
                );
                assert_eq!(
                    String::from_utf8(headers)?,
                    expected_headers,
                    "headers at {at}"
                );
                assert_eq!(
                    String::from_utf8(header_index)?,
                    expected_header_index,
                    "header .index at {at}"
                );
            }
        }

        Ok(())
    }
}
