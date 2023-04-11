use anyhow::{Context, Result};
use clap::ArgAction;
use clap::Parser;
use nale::align::bounded::structs::CloudSearchParams;
use nale::align::needleman_wunsch::{needleman_wunsch, SimpleTraceStep};
use nale::output::output_tabular::write_tabular_output;
use nale::pipelines::{pipeline_bounded, BoundedPipelineParams, DebugParams, Seed};
use nale::structs::hmm::parse_hmms_from_p7hmm_file;
use nale::structs::{Profile, Sequence};
use std::collections::HashMap;
use std::env;
use std::fs::{create_dir_all, File};
use std::io::{BufRead, BufReader, Read, Seek, SeekFrom};
use std::path::PathBuf;
use std::process::Command;
use std::time::Instant;
use thiserror::Error;

#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
pub struct Args {
    /// Query stockholm file
    query: String,
    /// Target fasta file
    target: String,
    /// Path for alignment output
    #[arg(long)]
    ali_out: Option<String>,
    /// Path for tabular output
    #[arg(long)]
    table_out: Option<String>,
    /// Write debugging information
    #[arg(long, action = ArgAction::SetTrue)]
    debug: Option<bool>,
    /// Allow output files to be overwritten
    #[arg(long, action = ArgAction::SetTrue)]
    allow_overwrite: Option<bool>,
}

impl Args {
    pub fn params_bounded(&self) -> BoundedPipelineParams {
        BoundedPipelineParams {
            allow_overwrite: match self.allow_overwrite {
                Some(val) => val,
                None => false,
            },
            cloud_search_params: CloudSearchParams::default(),
            debug_params: DebugParams {
                // TODO: parameterize this
                debug_path: PathBuf::from("./nale-debug"),
                write_matrices: match self.debug {
                    Some(val) => val,
                    None => false,
                },
                write_bounds: match self.debug {
                    Some(val) => val,
                    None => false,
                },
                write_trace: match self.debug {
                    Some(val) => val,
                    None => false,
                },
            },
        }
    }
}

#[derive(Error, Debug)]
#[error("no profile to profile map")]
struct ProfilesNotMappedError;

#[derive(Error, Debug)]
#[error("command exited without success")]
struct CommandExitStatusError;

trait CommandExt {
    fn run(&mut self) -> Result<()>;
}

impl CommandExt for Command {
    fn run(&mut self) -> Result<()> {
        let output = self.output().context("failed to start command")?;

        match output.status.success() {
            true => Ok(()),
            false => {
                let stdout = std::str::from_utf8(&output.stdout)
                    .context("failed to convert sdtout to UTF8")?;
                let stderr = std::str::from_utf8(&output.stderr)
                    .context("failed to convert sdterr to UTF8")?;
                println!("stdout: {stdout}");
                println!("stderr: {stderr}");
                Err(CommandExitStatusError.into())
            }
        }
    }
}

fn run_mmseqs(args: &Args, mmseqs_file_paths: &MmseqsFilePaths) -> Result<()> {
    Command::new("mmseqs")
        .arg("convertmsa")
        .arg(&args.query)
        .arg(&mmseqs_file_paths.query_msa_db)
        .run()?;

    Command::new("mmseqs")
        .arg("msa2profile")
        .arg(&mmseqs_file_paths.query_msa_db)
        .arg(&mmseqs_file_paths.query_db)
        .args(["--match-mode", "1"])
        .run()?;

    Command::new("mmseqs")
        .arg("createdb")
        .arg(&args.target)
        .arg(&mmseqs_file_paths.target_db)
        .run()?;

    println!("mmseqs prefilter");

    Command::new("mmseqs")
        .arg("prefilter")
        .arg(&mmseqs_file_paths.query_db)
        .arg(&mmseqs_file_paths.target_db)
        .arg(&mmseqs_file_paths.prefilter_db)
        // TODO: parameterize threads
        // .args(["--threads", "1"])
        // -k INT                    k-mer length (0: automatically set to optimum) [0]
        // TODO: I don't think we want to use this
        // .args(["-k", "7"])
        // --k-score INT             k-mer threshold for generating similar k-mer lists [2147483647]
        .args(["--k-score", "80"])
        // --min-ungapped-score INT  Accept only matches with ungapped alignment score above threshold [15]
        .args(["--min-ungapped-score", "15"])
        // --max-seqs INT            Maximum results per query sequence allowed to pass the prefilter (affects sensitivity) [300]
        .args(["--max-seqs", "1000"])
        .run()?;

    println!("mmseqs align");

    Command::new("mmseqs")
        .arg("align")
        .arg(&mmseqs_file_paths.query_db)
        .arg(&mmseqs_file_paths.target_db)
        .arg(&mmseqs_file_paths.prefilter_db)
        .arg(&mmseqs_file_paths.align_db)
        // TODO: parameterize threads
        // .args(["--threads", "1"])
        // -e DOUBLE      List matches below this E-value (range 0.0-inf) [1.000E-03]
        .args(["-e", "1e-2"])
        // --alt-ali INT  Show up to this many alternative alignments [0]
        .args(["--alt-ali", "0"])
        .args(["-a", "1"])
        .run()?;

    Command::new("mmseqs")
        .arg("convertalis")
        .arg(&mmseqs_file_paths.query_db)
        .arg(&mmseqs_file_paths.target_db)
        .arg(&mmseqs_file_paths.align_db)
        .arg(&mmseqs_file_paths.results)
        .args([
            "--format-output",
            "query,target,qstart,qend,tstart,tend,evalue",
        ])
        .run()?;

    Ok(())
}

fn extract_mmseqs_profile_consensus_sequence(
    mmseqs_file_paths: &MmseqsFilePaths,
) -> Result<HashMap<String, Sequence>> {
    let mut offsets_and_lengths: Vec<(usize, usize)> = vec![];
    let mut accession_numbers: Vec<String> = vec![];

    let query_db_h_index_file = File::open(&mmseqs_file_paths.query_db_h_index)
        .context("failed to open queryDB_h.index")?;

    let reader = BufReader::new(query_db_h_index_file);
    for line in reader.lines() {
        match line {
            Ok(line) => {
                let tokens: Vec<&str> = line.split_whitespace().collect();

                let offset = tokens[1].parse::<usize>()?;
                let length = tokens[2].parse::<usize>()?;
                offsets_and_lengths.push((offset, length));
            }
            Err(e) => {
                return Err(e).context("failed to parse line in queryDB_h.index");
            }
        }
    }

    let mut query_db_h_file =
        File::open(&mmseqs_file_paths.query_db_h).context("failed to open queryDB_h")?;

    for (offset, length) in &offsets_and_lengths {
        let mut buffer = vec![0; *length];
        query_db_h_file.seek(SeekFrom::Start(*offset as u64))?;
        query_db_h_file.read_exact(&mut buffer)?;

        let mut accession_string: Option<String> = None;
        for (buf_idx, byte) in buffer.iter().enumerate() {
            if byte.is_ascii_whitespace() {
                accession_string = Some(
                    std::str::from_utf8(&buffer[0..buf_idx])
                        .context("failed to create accession string")?
                        .to_string(),
                );
                break;
            }
        }

        match accession_string {
            Some(accession) => accession_numbers.push(accession),
            None => {
                panic!()
            }
        }
    }

    let query_db_index_file =
        File::open(&mmseqs_file_paths.query_db_index).context("failed to open queryDB.index")?;

    let reader = BufReader::new(query_db_index_file);
    for line in reader.lines() {
        match line {
            Ok(line) => {
                let tokens: Vec<&str> = line.split_whitespace().collect();

                let line_idx = tokens[0].parse::<usize>()?;
                let offset = tokens[1].parse::<usize>()?;
                let length = tokens[2].parse::<usize>()?;
                offsets_and_lengths[line_idx] = (offset, length);
            }
            Err(e) => {
                return Err(e).context("failed to parse line in queryDB.index");
            }
        }
    }

    let mut sequence_map: HashMap<String, Sequence> = HashMap::new();

    let mut query_db_file =
        File::open(&mmseqs_file_paths.query_db).context("failed to open queryDB")?;

    for (seq_idx, (offset, length)) in offsets_and_lengths.iter().enumerate() {
        let mut buffer = vec![0; *length];
        query_db_file.seek(SeekFrom::Start(*offset as u64))?;
        query_db_file.read_exact(&mut buffer)?;

        let mut consensus_digital_bytes: Vec<u8> = vec![];

        for byte_chunk in buffer.chunks(23) {
            if byte_chunk.len() == 23 {
                consensus_digital_bytes.push(byte_chunk[21]);
            }
        }

        sequence_map.insert(
            accession_numbers[seq_idx].clone(),
            Sequence::from_digital(&consensus_digital_bytes)?,
        );
    }

    Ok(sequence_map)
}

struct MmseqsFilePaths {
    query_msa_db: PathBuf,
    query_db: PathBuf,
    query_db_index: PathBuf,
    query_db_h: PathBuf,
    query_db_h_index: PathBuf,
    target_db: PathBuf,
    prefilter_db: PathBuf,
    align_db: PathBuf,
    results: PathBuf,
}

impl MmseqsFilePaths {
    fn new(root_path: &PathBuf) -> Self {
        Self {
            query_msa_db: root_path.join("msaDB"),
            query_db: root_path.join("queryDB"),
            query_db_index: root_path.join("queryDB.index"),
            query_db_h: root_path.join("queryDB_h"),
            query_db_h_index: root_path.join("queryDB_h.index"),
            target_db: root_path.join("targetDB"),
            prefilter_db: root_path.join("prefilterDB"),
            align_db: root_path.join("alignDB"),
            results: root_path.join("results.tsv"),
        }
    }
}

fn main() -> Result<()> {
    let args = Args::parse();

    let root_path = env::current_dir().unwrap().join("tmp");
    create_dir_all(&root_path)?;

    let mmseqs_file_paths = MmseqsFilePaths::new(&root_path);

    run_mmseqs(&args, &mmseqs_file_paths)?;
    let mmseqs_consensus_map = extract_mmseqs_profile_consensus_sequence(&mmseqs_file_paths)?;

    let mmseqs_results_file = File::open(&mmseqs_file_paths.results)?;
    let mmseqs_results_buf_reader = BufReader::new(mmseqs_results_file);

    // println!("hmmbuild...");
    let hmm_path = root_path.join("query.hmm");
    Command::new("hmmbuild")
        .arg(&hmm_path)
        .arg(&args.query)
        .run()?;

    println!("loading phmms...");

    let hmms = parse_hmms_from_p7hmm_file(hmm_path.to_str().unwrap())?;
    // let hmms = parse_hmms_from_p7hmm_file("./resources/pfam/pfam.hmm")?;

    let mut accession_to_name: HashMap<String, String> = HashMap::new();
    let mut name_to_accession: HashMap<String, String> = HashMap::new();
    for hmm in &hmms {
        accession_to_name.insert(hmm.header.accession_number.clone(), hmm.header.name.clone());
        name_to_accession.insert(hmm.header.name.clone(), hmm.header.accession_number.clone());
    }

    let profiles: Vec<Profile> = hmms.iter().map(|hmm| Profile::new(hmm)).collect();

    let mut profile_idx_maps_by_accession: HashMap<String, Vec<usize>> = HashMap::new();

    for profile in &profiles {
        let accession = name_to_accession.get(&profile.name).unwrap();
        let mmseqs_consensus = mmseqs_consensus_map.get(accession).unwrap();
        let p7_consensus = Sequence::from_utf8(&profile.consensus_sequence[1..])?;
        let trace = needleman_wunsch(mmseqs_consensus, &p7_consensus);

        let mut mmseqs_to_p7: Vec<usize> = vec![0; mmseqs_consensus.length + 1];

        let mut mmseqs_idx: usize = 0;
        let mut p7_idx: usize = 0;
        for step in &trace {
            match step {
                SimpleTraceStep::Diagonal => {
                    mmseqs_idx += 1;
                    p7_idx += 1;
                }
                SimpleTraceStep::Up => {
                    mmseqs_idx += 1;
                }
                SimpleTraceStep::Left => {
                    p7_idx += 1;
                }
            }
            mmseqs_to_p7[mmseqs_idx] = p7_idx;
        }
        debug_assert_eq!(mmseqs_idx, mmseqs_consensus.length);
        debug_assert_eq!(p7_idx, p7_consensus.length);

        profile_idx_maps_by_accession.insert(accession.clone(), mmseqs_to_p7);
    }

    let mut profile_seeds: HashMap<String, Vec<Seed>> = HashMap::new();

    for line in mmseqs_results_buf_reader.lines() {
        if let Ok(line) = line {
            let tokens: Vec<&str> = line.split_whitespace().collect();
            let accession = tokens[0];

            let profile_name = accession_to_name.get(accession).unwrap();

            let seeds = match profile_seeds.get_mut(profile_name) {
                Some(seeds) => seeds,
                None => {
                    profile_seeds.insert(profile_name.clone(), vec![]);
                    profile_seeds.get_mut(profile_name).unwrap()
                }
            };

            let profile_idx_map = profile_idx_maps_by_accession
                .get(accession)
                .ok_or(ProfilesNotMappedError)?;

            let target_name = tokens[1].to_string();
            let target_start = tokens[4].parse::<usize>()?;
            let target_end = tokens[5].parse::<usize>()?;
            let profile_start = tokens[2].parse::<usize>()?;
            let profile_end = tokens[3].parse::<usize>()?;

            // println!("{}, {}", target_start, profile_start);
            // println!("{}, {}", target_end, profile_end);

            seeds.push(Seed {
                target_name,
                target_start,
                target_end,
                profile_start: profile_idx_map[profile_start].max(1),
                profile_end: profile_idx_map[profile_end],
            })
        }
    }

    let mut profile_map: HashMap<String, Profile> = HashMap::new();

    for profile in profiles {
        profile_map.insert(profile.name.clone(), profile);
    }

    let targets = Sequence::amino_from_fasta(&args.target)?;
    let mut target_map: HashMap<String, Sequence> = HashMap::new();
    for target in targets {
        target_map.insert(target.name.clone(), target);
    }

    let params_bounded = args.params_bounded();
    let now = Instant::now();

    println!("nale start");
    let alignments_bounded = pipeline_bounded(
        &mut profile_map,
        &mut target_map,
        &profile_seeds,
        &params_bounded,
    )?;

    println!("{}µs", now.elapsed().as_micros());

    write_tabular_output(&alignments_bounded, &mut File::create("./results.out")?)?;
    Ok(())
}
