use anyhow::{Context, Result};
use clap::ArgAction;
use clap::Parser;
use nale::align::bounded::structs::CloudSearchParams;
use nale::output::output_tabular::write_tabular_output;
use nale::pipelines::{
    pipeline_bounded, pipeline_naive, BoundedPipelineParams, NaivePipelineParams, Seed,
};
use nale::structs::hmm::parse_hmms_from_p7hmm_file;
use nale::structs::{Profile, Sequence};
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::PathBuf;
use std::process::{Command, ExitStatus};
use std::time::Instant;
use std::{env, io};
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
            write_debug: match self.debug {
                Some(val) => val,
                None => false,
            },
            allow_overwrite: match self.allow_overwrite {
                Some(val) => val,
                None => false,
            },
            // TODO: parameterize this
            root_debug_dir_path: PathBuf::from("./nale-debug"),
            cloud_search_params: CloudSearchParams::default(),
        }
    }
}

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

fn main() -> Result<()> {
    let args = Args::parse();
    let now = Instant::now();

    let mut temp_dir_path = env::current_dir()?.join("tmp");

    let mmseqs_query_msa_db_path = temp_dir_path.join("msaDB");
    let mmseqs_query_db_path = temp_dir_path.join("queryDB");
    let mmseqs_target_db_path = temp_dir_path.join("targetDB");
    let mmseqs_prefilter_db_path = temp_dir_path.join("prefilterDB");
    let mmseqs_align_db_path = temp_dir_path.join("alignDB");
    let mmseqs_results_path = temp_dir_path.join("results.tsv");

    // Command::new("mmseqs")
    //     .arg("convertmsa")
    //     .arg(&args.query)
    //     .arg(&mmseqs_query_msa_db_path)
    //     .run()?;
    //
    // Command::new("mmseqs")
    //     .arg("msa2profile")
    //     .arg(&mmseqs_query_msa_db_path)
    //     .arg(&mmseqs_query_db_path)
    //     .args(["--match-mode", "1"])
    //     .run()?;
    //
    // Command::new("mmseqs")
    //     .arg("createdb")
    //     .arg(&args.target)
    //     .arg(&mmseqs_target_db_path)
    //     .run()?;
    //
    // println!("mmseqs prefilter");
    //
    // Command::new("mmseqs")
    //     .arg("prefilter")
    //     .arg(&mmseqs_query_db_path)
    //     .arg(&mmseqs_target_db_path)
    //     .arg(&mmseqs_prefilter_db_path)
    //     // .args(["-v", "1"])
    //     // .args(["--threads", "1"])
    //     // -k INT                      k-mer length (0: automatically set to optimum) [0]
    //     // .args(["-k", "7"])
    //     // --k-score INT               k-mer threshold for generating similar k-mer lists [2147483647]
    //     // .args(["--k-score", "80"])
    //     // --min-ungapped-score INT    Accept only matches with ungapped alignment score above threshold [15]
    //     // .args(["--min-ungapped-score", "15"])
    //     // --max-seqs INT              Maximum results per query sequence allowed to pass the prefilter (affects sensitivity) [300]
    //     // .args(["--max-seqs", "1000"])
    //     .run()?;
    //
    // println!("mmseqs align");
    //
    // Command::new("mmseqs")
    //     .arg("align")
    //     .arg(&mmseqs_query_db_path)
    //     .arg(&mmseqs_target_db_path)
    //     .arg(&mmseqs_prefilter_db_path)
    //     .arg(&mmseqs_align_db_path)
    //     // .args(["-v", "1"])
    //     // .args(["--threads", "1"])
    //     // -e DOUBLE                     List matches below this E-value (range 0.0-inf) [1.000E-03]
    //     .args(["-e", "1e-2"])
    //     // --alt-ali INT                 Show up to this many alternative alignments [0]
    //     .args(["--alt-ali", "0"])
    //     // -a BOOL                       Add backtrace string (convert to alignments with mmseqs convertalis module) [0]
    //     .args(["-a", "1"])
    //     .run()?;
    //
    // Command::new("mmseqs")
    //     .arg("convertalis")
    //     .arg(&mmseqs_query_db_path)
    //     .arg(&mmseqs_target_db_path)
    //     .arg(&mmseqs_align_db_path)
    //     .arg(&mmseqs_results_path)
    //     .args([
    //         "--format-output",
    //         "query,target,qstart,qend,tstart,tend,evalue",
    //     ])
    //     .run()?;

    let mmseqs_results_file = File::open(&mmseqs_results_path)?;
    let mmseqs_results_buf_reader = BufReader::new(mmseqs_results_file);

    // let hmm_path = temp_dir_path.join("query.hmm");
    // Command::new("hmmbuild")
    //     .arg(&hmm_path)
    //     .arg(&args.query)
    //     .run()?;
    //
    // println!("{}µs", now.elapsed().as_micros());

    println!("load phmm");
    // let hmms = parse_hmms_from_p7hmm_file(hmm_path.to_str().unwrap())?;
    let hmms = parse_hmms_from_p7hmm_file("./resources/pfam/pfam.hmm")?;

    let mut accession_to_name: HashMap<String, String> = HashMap::new();
    for hmm in &hmms {
        accession_to_name.insert(hmm.header.accession_number.clone(), hmm.header.name.clone());
    }

    let mut profile_seeds: HashMap<String, Vec<Seed>> = HashMap::new();
    for line in mmseqs_results_buf_reader.lines() {
        if let Ok(line) = line {
            let tokens: Vec<&str> = line.split_whitespace().collect();
            let profile_name = accession_to_name.get(tokens[0]).unwrap();
            let seeds = match profile_seeds.get_mut(profile_name) {
                Some(seeds) => seeds,
                None => {
                    profile_seeds.insert(profile_name.clone(), vec![]);
                    profile_seeds.get_mut(profile_name).unwrap()
                }
            };
            seeds.push(Seed {
                target_name: tokens[1].to_string(),
                target_start: tokens[4].parse::<usize>()?,
                target_end: tokens[5].parse::<usize>()?,
                profile_start: tokens[2].parse::<usize>()?,
                profile_end: tokens[3].parse::<usize>()?,
            })
        }
    }

    let profiles: Vec<Profile> = hmms.iter().map(|hmm| Profile::new(hmm)).collect();
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
