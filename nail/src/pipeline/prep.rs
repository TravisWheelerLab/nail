use crate::args::{guess_query_format_from_query_file, FileFormat};
use crate::cli::CommonArgs;
use crate::extension_traits::CommandExt;
use crate::pipeline::InvalidFileFormatError;

use std::fs::{self, create_dir_all};
use std::io::{BufWriter, Read, Write};
use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::{Arc, Mutex};

use anyhow::{Context, Result};
use clap::Args;
use libnail::structs::Sequence;
use rayon::iter::{IndexedParallelIterator, IntoParallelIterator, ParallelIterator};

#[derive(Args, Debug, Clone, Default)]
pub struct PrepDirArgs {
    /// Where MMseqs2 intermediate files are placed
    #[arg(
        short = 'p',
        long = "prep",
        default_value = "./prep/",
        value_name = "PATH"
    )]
    pub path: PathBuf,
}

impl PrepDirArgs {
    /// Produce a path to the query P7 HMM file
    pub fn prep_query_hmm_path(&self) -> PathBuf {
        self.path.join("query.hmm")
    }
    /// Produce a path to the MMseqs2 MSA database.
    ///
    /// This is created if a stockholm file query is provided.
    pub fn mmseqs_msa_db_path(&self) -> PathBuf {
        self.path.join("msaDB")
    }
    /// Produce a path to the MMseqs2 query database.
    ///
    /// If a fasta target was provided, this will be a sequence database.
    /// If a stockholm target was provided, this will be a profile database.
    pub fn mmseqs_query_db_path(&self) -> PathBuf {
        self.path.join("queryDB")
    }
    /// Produce a path to the MMseqs2 query database dbtype file.
    ///
    /// This file holds a byte that describes the original query file format.
    pub fn mmseqs_query_dbtype_path(&self) -> PathBuf {
        self.path.join("queryDB.dbtype")
    }
    /// Produce a path to the MMseqs2 query database index
    pub fn mmseqs_query_db_index_path(&self) -> PathBuf {
        self.path.join("queryDB.index")
    }
    /// Produce a path to the MMseqs2 query database h file
    pub fn mmseqs_query_db_h_path(&self) -> PathBuf {
        self.path.join("queryDB_h")
    }
    /// Produce a path to the MMseqs2 query database h file index
    pub fn mmseqs_query_db_h_index_path(&self) -> PathBuf {
        self.path.join("queryDB_h.index")
    }
    /// Produce a path to the MMseqs2 target database.
    ///
    /// This will always be a sequence database.
    pub fn mmseqs_target_db_path(&self) -> PathBuf {
        self.path.join("targetDB")
    }
    /// Produce a path to the MMseqs2 prefilter database.
    ///
    /// This is the result of running `mmseqs prefilter` on the query and target databases.
    pub fn mmseqs_prefilter_db_path(&self) -> PathBuf {
        self.path.join("prefilterDB")
    }
    /// Produce a path to the MMseqs2 alignment database.
    ///
    /// This is the result of running `mmseqs align` on the query, target, and prefilter databases.
    pub fn mmseqs_align_db_path(&self) -> PathBuf {
        self.path.join("alignDB")
    }
    /// Produce a path to the MMseqs2 alignment output.
    ///
    /// This is the result of running `mmseqs convertalis` on the query, target, and align databases.
    pub fn mmseqs_align_tsv_path(&self) -> PathBuf {
        self.path.join("align.tsv")
    }
}

#[derive(Args, Debug, Clone)]
pub struct PrepArgs {
    /// Query file
    #[arg(value_name = "QUERY.[fasta:sto]")]
    pub query_path: PathBuf,
    /// Target file
    #[arg(value_name = "TARGET.fasta")]
    pub target_path: PathBuf,
    /// Don't build a profile HMM with the input MSA
    #[arg(long, action)]
    pub skip_hmmbuild: bool,

    /// Provides paths to the prep dir and the files placed in it
    #[command(flatten)]
    pub prep_dir: PrepDirArgs,

    /// Arguments that are common across all nail subcommands
    #[command(flatten)]
    pub common_args: CommonArgs,
}

pub fn prep(args: &PrepArgs) -> Result<()> {
    let query_format = guess_query_format_from_query_file(&args.query_path)?;
    create_dir_all(&args.prep_dir.path)?;

    match query_format {
        FileFormat::Fasta => {
            Command::new("mmseqs")
                .arg("createdb")
                .arg(&args.query_path)
                .arg(&args.prep_dir.mmseqs_query_db_path())
                .run()?;

            if !args.skip_hmmbuild {
                build_hmm_from_fasta(
                    &args.query_path,
                    &args.prep_dir.path.join("tmp.fasta"),
                    &args.prep_dir.prep_query_hmm_path(),
                    &args.prep_dir.path.join("tmp.hmm"),
                    args.common_args.num_threads,
                )?;
            }
        }
        FileFormat::Stockholm => {
            Command::new("mmseqs")
                .arg("convertmsa")
                // this flag should force the profile to be labeled
                // using the actual name rather than the accession number
                .args(["--identifier-field", "0"])
                .arg(&args.query_path)
                .arg(&args.prep_dir.mmseqs_msa_db_path())
                .run()?;

            Command::new("mmseqs")
                .arg("msa2profile")
                .arg(&args.prep_dir.mmseqs_msa_db_path())
                .arg(&args.prep_dir.mmseqs_query_db_path())
                .args(["--threads", &args.common_args.num_threads.to_string()])
                // --match-mode INT       0: Columns that have a residue in the first sequence are kept,
                //                        1: columns that have a residue in --match-ratio of all sequences
                //                           are kept [0]
                .args(["--match-mode", "1"])
                .run()?;

            if !args.skip_hmmbuild {
                build_hmm_from_stockholm(
                    &args.query_path,
                    &args.prep_dir.prep_query_hmm_path(),
                    args.common_args.num_threads,
                )?;
            }
        }
        ref format => {
            return Err(InvalidFileFormatError {
                format: format.clone(),
            })
            .context("invalid query file format in nail prep")
        }
    }

    Command::new("mmseqs")
        .arg("createdb")
        .arg(&args.target_path)
        .arg(&args.prep_dir.mmseqs_target_db_path())
        .run()?;

    Ok(())
}

pub fn build_hmm_from_stockholm(
    stockholm_path: &impl AsRef<Path>,
    hmm_path: &impl AsRef<Path>,
    num_threads: usize,
) -> Result<()> {
    Command::new("hmmbuild")
        .args(["--cpu", &num_threads.to_string()])
        .arg(hmm_path.as_ref())
        .arg(stockholm_path.as_ref())
        .run()?;

    Ok(())
}

pub fn build_hmm_from_fasta(
    query_fasta_path: &impl AsRef<Path>,
    temp_fasta_path: &impl AsRef<Path>,
    query_hmm_path: &impl AsRef<Path>,
    temp_hmm_path: &impl AsRef<Path>,
    num_threads: usize,
) -> anyhow::Result<()> {
    let query_fasta_path = query_fasta_path.as_ref();
    let query_hmm_path = query_hmm_path.as_ref();

    let query_seqs = Sequence::amino_from_fasta(query_fasta_path).with_context(|| {
        format!(
            "failed to parse query fasta: {}",
            query_fasta_path.to_string_lossy()
        )
    })?;

    // this is the query hmm file we are going to build
    let query_hmm_file = std::fs::File::create(query_hmm_path)?;
    let hmm_buf_writer = BufWriter::new(query_hmm_file);

    fs::create_dir(temp_hmm_path)?;
    fs::create_dir(temp_fasta_path)?;

    let pb = indicatif::ProgressBar::new(query_seqs.len() as u64);
    pb.set_style(indicatif::ProgressStyle::with_template("{spinner:.green} [{elapsed_precise}] [{bar:.cyan/blue}] {human_pos}/{human_len} ({eta})")
        .unwrap()
        .progress_chars("#>-"));

    #[derive(Clone)]
    struct ThreadData {
        progress_bar: Arc<Mutex<indicatif::ProgressBar>>,
        hmm_writer: Arc<Mutex<BufWriter<std::fs::File>>>,
        hmm_bytes: Vec<u8>,
        hmm_path: PathBuf,
        fasta_path: PathBuf,
    }

    let init = ThreadData {
        progress_bar: Arc::new(Mutex::new(pb)),
        hmm_writer: Arc::new(Mutex::new(hmm_buf_writer)),
        hmm_bytes: vec![],
        hmm_path: temp_hmm_path.as_ref().to_path_buf(),
        fasta_path: temp_fasta_path.as_ref().to_path_buf(),
    };

    query_seqs
        .into_par_iter()
        .enumerate()
        .try_for_each_with(init, |data, (idx, seq)| {
            data.hmm_bytes.clear();

            let temp_fasta_path = data.fasta_path.join(idx.to_string());
            let temp_hmm_path = data.hmm_path.join(idx.to_string());

            let mut temp_fasta_file = std::fs::File::create(&temp_fasta_path)
                .context("failed to open temporary fasta file")?;
            writeln!(temp_fasta_file, "{}", seq)
                .context("failed to write to temporary fasta file")?;

            Command::new("hmmbuild")
                .args(["-n", &seq.name])
                .arg("--amino")
                .arg(&temp_hmm_path)
                .arg(&temp_fasta_path)
                .run()
                .context("failed to run hmmbuild")?;

            let mut temp_hmm_file =
                std::fs::File::open(temp_hmm_path).context("failed to open temporary hmm file")?;

            temp_hmm_file
                .read_to_end(&mut data.hmm_bytes)
                .context("failed to read temporary hmm file")?;

            let pb = data.progress_bar.lock().unwrap();
            pb.inc(1);

            let mut writer = data.hmm_writer.lock().unwrap();
            writer
                .write_all(&data.hmm_bytes)
                .context("failed to write to query hmm file")
        })
        .context("failed to build hmm from fasta")?;

    Ok(())
}
