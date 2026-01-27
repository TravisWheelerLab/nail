use self::consts::*;

use std::{
    fs::File,
    io::{BufRead, BufReader, Read, Seek, Write},
    path::{Path, PathBuf},
    process::Command,
};

use libnail::{align::Nats, alphabet::UTF8_TO_DIGITAL_AMINO, structs::Profile};

use crate::{
    args::SearchArgs,
    io::{Database, Fasta, ReadSeekExt, ReadState},
    util::{CommandExt, PathBufExt},
};

use anyhow::Context;
use anyhow::{anyhow, bail};
use regex::Regex;

pub mod consts {
    pub const AMINO_DBTYPE: &[u8] = &[0, 0, 0, 0];
    pub const PROFILE_DBTYPE: &[u8] = &[2, 0, 0, 0];
    pub const ALIGN_DBTYPE: &[u8] = &[5, 0, 0, 0];
    pub const PREFILTER_DBTYPE: &[u8] = &[7, 0, 0, 0];
    pub const GENERIC_DBTYPE: &[u8] = &[12, 0, 0, 0];
    #[cfg(not(target_os = "windows"))]
    pub const BLOSUM_62: &str = include_str!("../mat/blosum62.mat");
    #[cfg(not(target_os = "windows"))]
    pub const BLOSUM_80: &str = include_str!("../mat/blosum80.mat");
    #[cfg(target_os = "windows")]
    pub const BLOSUM_62: &str = include_str!("..\\mat\\blosum62.mat");
    #[cfg(target_os = "windows")]
    pub const BLOSUM_80: &str = include_str!("..\\mat\\blosum80.mat");
}

pub enum ByteBuffer<'a> {
    Complete(&'a [u8]),
    Partial(&'a [u8], usize),
    Empty,
}

#[derive(Clone, Copy)]
struct Descriptor {
    file_idx: usize,
    offset: u64,
    length: u64,
}

pub struct PrefilterDb {
    file: File,
    paths: Vec<PathBuf>,
    descriptors: Vec<Descriptor>,
    checkpoints: Vec<u64>,
    buf: Vec<u8>,
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

        let mut descriptors = offsets
            .into_iter()
            .zip(lengths)
            .map(|(offset, length)| {
                let file_idx = match prefix_sum.binary_search(&offset) {
                    // this means the binary search found the exact offset
                    Ok(idx) => idx,
                    // this means the binary search landed between offsets
                    Err(idx) => idx.saturating_sub(1),
                };

                let relative_offset = offset - prefix_sum[file_idx];
                Descriptor {
                    file_idx,
                    offset: relative_offset,
                    length,
                }
            })
            .collect::<Vec<_>>();

        for desc in descriptors.iter_mut() {
            let mut file = File::open(&paths[desc.file_idx])?;
            file.seek(std::io::SeekFrom::Start(desc.offset + desc.length))?;

            let mut b = [0u8; 1];
            file.read_exact(&mut b)?;
            let byte = b[0];
            assert!(byte == 0);
        }

        let file = File::open(&paths[0])?;
        let checkpoints = vec![0; descriptors.len()];
        let buf = Vec::with_capacity(2 << 11);

        Ok(Self {
            file,
            paths,
            descriptors,
            checkpoints,
            buf,
        })
    }

    pub fn open_file(&mut self, file_idx: usize) -> anyhow::Result<()> {
        self.file = File::open(&self.paths[file_idx])?;
        Ok(())
    }

    pub fn get(&mut self, prf_idx: usize) -> anyhow::Result<&[u8]> {
        let desc = self.descriptors[prf_idx];

        if desc.length == 0 {
            return Ok(&[]);
        }

        self.open_file(desc.file_idx)?;

        self.file.seek(std::io::SeekFrom::Start(desc.offset))?;
        let mut taken = (&mut self.file).take(desc.length);

        self.buf.clear();
        self.buf.resize(desc.length as usize, 0);

        match taken.read_with_state(&mut self.buf)? {
            ReadState::Final(_) => Ok(&self.buf),
            _ => unreachable!(),
        }
    }

    pub fn next_n(&'_ mut self, prf_idx: usize, n: usize) -> anyhow::Result<ByteBuffer<'_>> {
        if n == 0 {
            return Ok(ByteBuffer::Empty);
        }

        let desc = self.descriptors[prf_idx];
        self.open_file(desc.file_idx)?;

        let checkpoint = &mut self.checkpoints[prf_idx];

        let start = desc.offset + *checkpoint;
        let end = desc.offset + desc.length;
        let remaining = end.saturating_sub(start);

        if remaining == 0 {
            return Ok(ByteBuffer::Empty);
        }

        self.file.seek(std::io::SeekFrom::Start(start))?;
        let mut taken = (&mut self.file).take(remaining);

        self.buf.clear();

        const CHUNK_SZ: usize = 2 << 12;
        let mut record_cnt = 0;
        loop {
            let old_len = self.buf.len();
            self.buf.resize(old_len + CHUNK_SZ, 0);

            let n_read = match taken.read_with_state(&mut self.buf[old_len..])? {
                ReadState::Reading(n) | ReadState::Final(n) => n,
                ReadState::Done => {
                    self.buf.truncate(old_len);
                    *checkpoint = desc.length;
                    return Ok(ByteBuffer::Partial(&self.buf, record_cnt));
                }
            };

            self.buf.truncate(old_len + n_read);

            for (pos, &b) in self.buf.iter().enumerate().skip(old_len) {
                if b == b'\n' {
                    record_cnt += 1;
                }

                if record_cnt == n {
                    *checkpoint += (pos + 1) as u64;
                    self.buf.truncate(pos + 1);
                    return Ok(ByteBuffer::Complete(&self.buf));
                }
            }
        }
    }
}

#[allow(dead_code)]
#[derive(Clone, Copy)]
pub enum MmseqsScoreModel {
    Profile,
    Blosum62,
    Blosum80,
}

impl MmseqsScoreModel {
    pub fn write(&self, dir: impl AsRef<Path>) -> anyhow::Result<Option<PathBuf>> {
        let dir = dir.as_ref();

        if !dir.is_dir() {
            bail!("path: {} is not a directory", dir.to_string_lossy());
        }

        let (mat_str, mat_path) = match self {
            MmseqsScoreModel::Profile => return Ok(None),
            MmseqsScoreModel::Blosum62 => (BLOSUM_62, dir.join("blosum62.out")),
            MmseqsScoreModel::Blosum80 => (BLOSUM_80, dir.join("blosum80.out")),
        };

        mat_path.open(true)?.write_all(mat_str.as_bytes())?;

        Ok(Some(mat_path))
    }
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
        let dir = dir.as_ref().to_path_buf();
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

    let dir = db_path.parent().unwrap();
    std::fs::create_dir_all(dir)?;

    db_dbtype_path.open(true)?.write_all(AMINO_DBTYPE)?;
    header_dbtype_path.open(true)?.write_all(GENERIC_DBTYPE)?;

    let mut db = db_path.open(true)?;
    let mut db_index = db_index_path.open(true)?;
    let mut db_header = header_path.open(true)?;
    let mut db_header_index = header_index_path.open(true)?;

    let mut db_offset = 0usize;
    let mut header_offset = 0usize;

    for (seq_count, seq) in sequences.values().enumerate() {
        db.write_all(&seq.utf8_bytes[1..])?;
        db.write_all(&[10u8, 0u8])?;

        let db_byte_length = seq.length + 2;
        writeln!(db_index, "{}\t{}\t{}", seq_count, db_offset, db_byte_length,)?;

        writeln!(db_header, "{}", seq.name)?;
        db_header.write_all(&[0u8])?;

        let header_byte_length = seq.name.len() + 2;
        writeln!(
            db_header_index,
            "{}\t{}\t{}",
            seq_count, header_offset, header_byte_length
        )?;

        db_offset += db_byte_length;
        header_offset += header_byte_length;
    }

    Ok(())
}

pub fn write_mmseqs_profile_database(
    profiles: impl Iterator<Item = Profile>,
    path: impl AsRef<Path>,
) -> anyhow::Result<()> {
    let db_path = path.as_ref().to_owned();
    let db_name = db_path.file_name().unwrap().to_str().unwrap();
    let db_index_path = db_path.with_file_name(format!("{db_name}.index"));
    let db_dbtype_path = db_path.with_file_name(format!("{db_name}.dbtype"));

    let header_path = db_path.with_file_name(format!("{db_name}_h"));
    let header_index_path = db_path.with_file_name(format!("{db_name}_h.index"));
    let header_dbtype_path = db_path.with_file_name(format!("{db_name}_h.dbtype"));

    let dir = db_path.parent().unwrap();
    std::fs::create_dir_all(dir)?;

    db_dbtype_path.open(true)?.write_all(PROFILE_DBTYPE)?;
    header_dbtype_path.open(true)?.write_all(GENERIC_DBTYPE)?;

    let mut db = db_path.open(true)?;
    let mut db_index = db_index_path.open(true)?;
    let mut db_header = header_path.open(true)?;
    let mut db_header_index = header_index_path.open(true)?;

    let mut db_offset = 0usize;
    let mut header_offset = 0usize;

    for (prf_cnt, prf) in profiles.enumerate() {
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
        writeln!(db_header, "{}", prf.name)?;
        db_header.write_all(&[0u8])?;

        // +1 for the 0 byte, +1 for the newline
        let header_byte_length = prf.name.len() + 2;

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
    let mut prefilter = Command::new("mmseqs");

    let qdb = query_db_path.as_ref();
    let tdb = target_db_path.as_ref();
    let pdb = prefilter_db_path.as_ref();

    let pdb_dir = pdb
        .parent()
        .context("failed to produce mmseqs prefilter DB directory path")?;

    if !pdb_dir
        .try_exists()
        .with_context(|| format!("failed to check existence of: {pdb_dir:?}"))?
    {
        std::fs::create_dir(pdb_dir).context("failed to create mmseqs prefilter DB directory")?;
    }

    match pdb.parent() {
        Some(dir) => std::fs::create_dir_all(dir)?,
        None => bail!("failed to create mmseqs prefilter DB directory"),
    }

    prefilter
        .arg("prefilter")
        .arg(qdb)
        .arg(tdb)
        .arg(pdb)
        .args(["--threads", &args.num_threads.to_string()])
        .args(["-k", &args.mmseqs_args.k.to_string()])
        .args(["-s", &args.mmseqs_args.s.to_string()])
        .args(["--max-seqs", &args.mmseqs_args.max_seqs.to_string()]);

    if let Some(v) = args.mmseqs_args.comp_bias_corr {
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

    if !adb_dir
        .try_exists()
        .with_context(|| format!("failed to check existence of: {adb_dir:?}"))?
    {
        std::fs::create_dir(adb_dir).context("failed to create mmseqs align DB directory")?;
    }

    let mut align = Command::new("mmseqs");
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

    if let Some(v) = args.mmseqs_args.comp_bias_corr {
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

    Command::new("mmseqs")
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
