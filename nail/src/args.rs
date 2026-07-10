use std::{
    io::{IsTerminal, Write},
    path::PathBuf,
    str::FromStr,
};

use anyhow::{anyhow, bail, Context};
use clap::{Args, Parser, Subcommand, ValueEnum};
use regex::Regex;

use crate::util::{term::*, PathExt};

#[derive(Parser)]
#[command(
    version,
    name = "nail",
    about = "Using MMseqs2 to find rough alignment seeds, perform bounded profile HMM sequence alignment"
)]
pub struct NailCli {
    #[command(subcommand)]
    pub command: NailSubCommands,
}

#[allow(clippy::large_enum_variant)]
#[derive(Subcommand)]
pub enum NailSubCommands {
    #[command(about = "Run nail's protein search pipeline")]
    Search(SearchArgs),
    #[command(subcommand, hide = true)]
    Dev(DevSubCommands),
}

#[allow(clippy::large_enum_variant)]
#[derive(Subcommand)]
pub enum DevSubCommands {
    #[command()]
    Play(SearchArgs),
    Search(SearchArgs),
    Mx(SearchArgs),
}

#[derive(Debug, Args)]
pub struct SearchArgs {
    /// The query database file
    #[arg(value_name = "QUERY.[fasta:hmm]")]
    pub query_path: PathBuf,

    /// The target database file
    #[arg(value_name = "TARGET.fasta")]
    pub target_path: PathBuf,

    /// The number of threads that nail will use
    #[arg(short = 't', default_value_t = 8usize, value_name = "N")]
    pub num_threads: usize,

    /// Print out pipeline summary statistics
    #[arg(short = 's', action)]
    pub print_summary_stats: bool,

    /// Don't write any tabular results, write alignments to stdout
    #[arg(short = 'x', action)]
    pub ali_to_stdout: bool,

    /// The path to the MMseqs2 binary that nail will use for seeding
    #[arg(long, value_name = "/path/to/mmseqs", default_value = "mmseqs")]
    pub mmseqs_path: PathBuf,

    #[command(flatten)]
    #[clap(next_help_heading = "File I/O options")]
    pub io_args: IoArgs,

    #[command(flatten)]
    #[clap(next_help_heading = "Pipeline options")]
    pub pipeline_args: PipelineArgs,

    /// Arguments that are passed to MMseqs2
    #[command(flatten)]
    #[clap(next_help_heading = "Seeding options")]
    pub seed_args: SeedArgs,

    #[command(flatten)]
    #[clap(next_help_heading = "Expert options")]
    pub expert_args: ExpertArgs,

    #[command(flatten)]
    #[clap(next_help_heading = "Dev options")]
    pub dev_args: DevArgs,
}

impl SearchArgs {
    pub fn validate(&mut self) -> anyhow::Result<()> {
        // ----------------
        // output precision
        // ----------------

        if let Some(p1) = self.dev_args.bit_p {
            let p2 = *libnail::output::output_tabular::BIT_P.get_or_init(|| p1);
            if p1 != p2 {
                bail!("failed to set BIT_P")
            }
        }

        if let Some(p1) = self.dev_args.f32_p {
            let p2 = *libnail::output::output_tabular::F32_P.get_or_init(|| p1);
            if p1 != p2 {
                bail!("failed to set F32_P")
            }
        }

        if let Some(p1) = self.dev_args.f64_p {
            let p2 = *libnail::output::output_tabular::F64_P.get_or_init(|| p1);
            if p1 != p2 {
                bail!("failed to set F64_P")
            }
        }

        match self.io_args.tmp_dir_path {
            Some(ref dir) => {
                dir.create_dir()?;
            }
            None => {
                let root_dir = PathBuf::from("./tmp-nail");
                root_dir.create_dir()?;

                let now = std::time::Instant::now();
                const TIMEOUT: u64 = 1;
                let dir = loop {
                    if now.elapsed().as_secs() >= TIMEOUT {
                        break Err(anyhow!("timeout: {TIMEOUT}s on Unix timestamping"));
                    }

                    let unix_ns = std::time::SystemTime::now()
                        .duration_since(std::time::UNIX_EPOCH)
                        .context("failed to produce Unix epoch time")?
                        .as_nanos();

                    let dir = root_dir.join(unix_ns.to_string());

                    match std::fs::create_dir(&dir) {
                        Ok(_) => break Ok(dir),
                        Err(e) if e.kind() == std::io::ErrorKind::AlreadyExists => {
                            // this probably almost never matters, but yielding
                            // should help the case where many, many threads
                            // happen to collide on the directory namespace
                            std::thread::yield_now();
                            continue;
                        }
                        Err(e) => break Err(e.into()),
                    }
                }
                .with_context(|| {
                    format!(
                        "failed to create temp database directory under: {:?}",
                        root_dir,
                    )
                })?;

                self.io_args.tmp_dir_path = Some(dir);
            }
        }

        {
            // quickly make sure we can write to all of the results paths

            if let Some(path) = &self.io_args.tbl_results_path {
                path.check_open(self.io_args.allow_overwrite)?;
            }

            if let Some(path) = &self.io_args.ali_results_path {
                path.check_open(self.io_args.allow_overwrite)?;
            }

            if let Some(path) = &self.io_args.seeds_output_path {
                path.check_open(self.io_args.allow_overwrite)?;
            }

            if let Some(path) = &self.dev_args.stats_results_path {
                path.check_open(self.io_args.allow_overwrite)?;
            }
        }

        if self.ali_to_stdout {
            self.io_args.tbl_results_path = None
        }

        if self.pipeline_args.only_seed && self.io_args.seeds_output_path.is_none() {
            self.io_args.seeds_output_path = Some(PathBuf::from_str("./seeds.tsv")?);
        }

        match self.seed_args.seed_mode {
            SeedMode::Static => {
                if self.seed_args.prog_n.is_some() {
                    bail!("the argument '{YELLOW}--prog-n{RESET}' cannot be used without '{YELLOW}--prog-seed{RESET}'")
                }
                if self.seed_args.prog_f.is_some() {
                    bail!("the argument '{YELLOW}--prog-f{RESET}' cannot be used without '{YELLOW}--prog-seed{RESET}'")
                }
            }
            SeedMode::Prog => {
                if self.seed_args.prog_n.is_none() {
                    bail!("the argument '{YELLOW}--prog-n{RESET}' is unset")
                }
                if self.seed_args.prog_f.is_none() {
                    bail!("the argument '{YELLOW}--prog-f{RESET}' is unset")
                }
            }
        }

        Ok(())
    }

    pub fn write(&self, out: &mut impl Write) -> anyhow::Result<()> {
        writeln!(
            out,
            "target: {}",
            self.target_path.to_str().unwrap_or_default()
        )?;
        writeln!(
            out,
            "query:  {}",
            self.query_path.to_str().unwrap_or_default()
        )?;
        writeln!(out)?;

        writeln!(out, "pipeline arguments:")?;
        writeln!(out, " ├─ mmseqs -k: {}", self.seed_args.k)?;
        writeln!(out, " ├─ mmseqs -s: {:.}", self.seed_args.s)?;
        writeln!(out, " ├─ mmseqs --max-seqs: {:.}", self.seed_args.max_seqs)?;
        match self.seed_args.seed_mode {
            SeedMode::Static => {
                writeln!(out, " ├─ seed-mode: static",)?;
            }

            SeedMode::Prog => {
                writeln!(out, " ├─ seed-mode: prog",)?;

                writeln!(
                    out,
                    "   ├─ n: {}",
                    self.seed_args.prog_n.unwrap_or_default()
                )?;
                writeln!(
                    out,
                    "   └─ f: {}",
                    self.seed_args.prog_f.unwrap_or_default()
                )?;
            }
        }

        writeln!(out, " ├─ α: {}", self.pipeline_args.alpha)?;
        writeln!(out, " ├─ β: {}", self.pipeline_args.beta)?;
        writeln!(out, " ├─ γ: {}", self.pipeline_args.gamma)?;
        writeln!(out, " ├─ S: {}", self.pipeline_args.seed_pvalue_threshold)?;
        writeln!(out, " ├─ C: {}", self.pipeline_args.cloud_pvalue_threshold)?;
        writeln!(
            out,
            " ├─ F: {}",
            self.pipeline_args.forward_pvalue_threshold
        )?;
        writeln!(out, " ├─ E: {}", self.pipeline_args.e_value_threshold)?;
        writeln!(
            out,
            " └─ Z: {}",
            self.expert_args.target_database_size.unwrap_or_default()
        )?;

        Ok(())
    }
}

#[derive(Args, Debug, Clone, Default)]
pub struct IoArgs {
    /// The file where tabular output will be written
    #[arg(long = "tbl-out", default_value = "results.tbl", value_name = "PATH")]
    pub tbl_results_path: Option<PathBuf>,

    /// The file where alignment output will be written
    #[arg(long = "ali-out", default_value = None, value_name = "PATH")]
    pub ali_results_path: Option<PathBuf>,

    /// A file containing pre-computed alignment seeds
    #[arg(long = "seeds", value_name = "PATH")]
    pub seeds_input_path: Option<PathBuf>,

    /// The file where alignment seeds will be written
    #[arg(long = "seeds-out", default_value = None, value_name = "PATH")]
    pub seeds_output_path: Option<PathBuf>,

    /// The directory where intermediate files will be placed
    #[arg(long = "tmp-dir", value_name = "PATH")]
    pub tmp_dir_path: Option<PathBuf>,

    /// Allow nail to overwrite files
    #[arg(short = 'X', long = "allow-overwrite", default_value_t = false)]
    pub allow_overwrite: bool,
}

#[derive(Args, Debug, Clone, Default)]
pub struct PipelineArgs {
    /// (cloud search) local score pruning threshold
    #[arg(
        short = 'A',
        default_value_t = 10.0,
        value_name = "X",
        verbatim_doc_comment
    )]
    pub alpha: f32,

    /// (cloud search) global score pruning threshold
    #[arg(
        short = 'B',
        default_value_t = 16.0,
        value_name = "X",
        verbatim_doc_comment
    )]
    pub beta: f32,

    /// (cloud search) at minimum, compute N anti-diagonals
    #[arg(
        short = 'G',
        default_value_t = 5,
        value_name = "N",
        verbatim_doc_comment
    )]
    pub gamma: usize,

    /// (cloud search) allow at most N attempts to recover disjoint clouds
    #[arg(
        short = 'a',
        default_value_t = 5,
        value_name = "N",
        verbatim_doc_comment
    )]
    pub cloud_max_join_attempts: usize,

    /// (cloud search) when cloud search fails, scale cloud search parameters (α, ɓ, γ) by 1.0 + X
    #[arg(
        short = 'f',
        default_value_t = 0.5,
        value_name = "X",
        verbatim_doc_comment
    )]
    pub cloud_param_scale_factor: f32,

    /// (seeding filter) filter hits with P-value > X
    #[arg(
        short = 'S',
        default_value_t = 1e-4,
        value_name = "X",
        verbatim_doc_comment
    )]
    pub seed_pvalue_threshold: f64,

    /// (cloud search filter) filter hits with P-value > X
    #[arg(
        short = 'C',
        default_value_t = 1e-2,
        value_name = "X",
        verbatim_doc_comment
    )]
    pub cloud_pvalue_threshold: f64,

    /// (forward filter) filter hits with P-value > X
    #[arg(
        short = 'F',
        default_value_t = 1e-4,
        value_name = "X",
        verbatim_doc_comment
    )]
    pub forward_pvalue_threshold: f64,

    /// (reporting filter) filter hits with E-value > X
    #[arg(
        short = 'E',
        default_value_t = 10.0,
        value_name = "X",
        verbatim_doc_comment
    )]
    pub e_value_threshold: f64,

    /// Produce alignment seeds and terminate
    #[arg(long = "only-seed", action)]
    pub only_seed: bool,
}

#[derive(Args, Debug, Clone, Default)]
pub struct ExpertArgs {
    /// Override the number of comparisons used for E-value calculation
    #[arg(short = 'Z', value_name = "N")]
    pub target_database_size: Option<usize>,

    /// Don't compute sequence composition bias score correction
    #[arg(long = "no-null2", action)]
    pub no_null_two: bool,
}

#[derive(Args, Debug, Clone, Default)]
pub struct DevArgs {
    /// Where to place stats output
    #[arg(long, value_name = "PATH", hide = true)]
    pub stats_results_path: Option<PathBuf>,

    /// Compute the full DP matrices
    #[arg(long, action, hide = true)]
    pub full_dp: bool,

    /// Formatting precision for Bits
    #[arg(long, hide = true)]
    pub bit_p: Option<usize>,

    /// Formatting precision for f32
    #[arg(long, hide = true)]
    pub f32_p: Option<usize>,

    /// Formatting precision for f64
    #[arg(long, hide = true)]
    pub f64_p: Option<usize>,
}

#[derive(Args, Debug, Clone, Default)]
pub struct SeedArgs {
    /// (MMseqs2 prefilter) k-mer length (0: automatically set to optimum)
    #[arg(
        long = "mmseqs-k",
        default_value_t = 6usize,
        value_name = "N",
        verbatim_doc_comment
    )]
    pub k: usize,

    /// (MMseqs2 prefilter) Sensitivity: 1.0 faster; 4.0 fast; 7.5 sensitive
    #[arg(
        long = "mmseqs-s",
        default_value_t = 10.0,
        value_name = "X",
        verbatim_doc_comment
    )]
    pub s: f32,

    /// (MMseqs2 prefilter): Maximum results per query sequence allowed to pass the prefilter
    ///
    ///
    #[arg(
        long = "mmseqs-max-seqs",
        default_value_t = 2147483647usize,
        default_value_if("seed_mode", "static", "300"),
        default_value_if("seed_mode", "prog", "2147483647"),
        value_name = "N",
        verbatim_doc_comment
    )]
    pub max_seqs: usize,

    /// Impose a limit on the number of seeds that nail will align per query
    #[arg(long = "max-seeds", value_name = "N", verbatim_doc_comment)]
    pub max_seeds: Option<usize>,

    /// (MMseqs2 Parameter) Correct for locally biased amino acid composition (range 0-1)
    #[arg(
        long = "mmseqs-comp-bias-corr",
        value_name = "N",
        hide = true,
        verbatim_doc_comment
    )]
    pub comp_bias_corr: Option<usize>,

    /// How nail will produce alignment seeds
    ///
    ///   static|0: run the MMseqs2 prefilter, then run MMseqs2 align on the entire prefilter results
    ///     prog|1: run the MMseqs2 prefilter, then progressively run MMseqs2 align following
    ///             parameterization of --prog-n and prog-f
    #[arg(
        long,
        default_value = "prog",
        value_name = "MODE",
        verbatim_doc_comment
    )]
    pub seed_mode: SeedMode,

    /// (prog seeding) the initial number of mmseqs alignments per query
    #[arg(
        long,
        default_value = "200",
        default_value_if("seed_mode", "static", None),
        default_value_if("seed_mode", "prog", "200"),
        value_name = "N",
        verbatim_doc_comment
    )]
    pub prog_n: Option<usize>,

    /// (prog seeding) the fraction of hits required to continue progressive seeding
    #[arg(
        long,
        default_value = "0.01",
        default_value_if("seed_mode", "static", None),
        default_value_if("seed_mode", "prog", "0.01"),
        value_name = "X",
        verbatim_doc_comment
    )]
    pub prog_f: Option<f32>,
}

#[derive(Default, Clone, Copy, Debug, ValueEnum)]
pub enum SeedMode {
    #[value(alias = "0")]
    Static,
    #[default]
    #[value(alias = "1")]
    Prog,
}

pub fn handle_clap_error(e: clap::Error) -> ! {
    match e.kind() {
        clap::error::ErrorKind::DisplayHelp => {
            let args: Vec<String> = std::env::args().collect();

            let is_terminal = std::io::stdout().is_terminal();
            let is_short_help_zebra = args[2] == "-hz";

            if is_terminal && is_short_help_zebra {
                print_help_summary_zebra(&e);
                std::process::exit(e.exit_code());
            } else {
                e.exit();
            }
        }
        _ => e.exit(),
    }
}

pub fn print_help_summary_zebra(e: &clap::Error) {
    let ansi_string = format!("{}", e.render().ansi());

    // ANSI codes escaped for use in regex search
    const BOLD_RE: &str = r"\x1b\[1m";
    const UNDERLINE_RE: &str = r"\x1b\[4m";
    const RESET_RE: &str = r"\x1b\[0m";

    // this should find any ANSI formatted clap --help line that is a header
    // i.e. this will find "usage", "arguments," "options," and then anything
    // set up in the CLI using the "(next_help_heading macro = <HEADER>)"
    let header_re = Regex::new(&format!(r"(?m)^{BOLD_RE}{UNDERLINE_RE}.*:{RESET_RE}")).unwrap();

    /// Given a regex and string, do a split-on-regex-inclusive search.
    ///
    /// returns: (splits, matches)
    fn split_and_match(re: &Regex, str: &str) -> (Vec<String>, Vec<String>) {
        let mut other = Vec::new();
        let mut matches = Vec::new();
        let mut last = 0;

        for m in re.find_iter(str) {
            other.push(str[last..m.start()].to_string());
            matches.push(m.as_str().to_string());
            last = m.end();
        }

        other.push(str[last..].to_string());

        (other, matches)
    }

    let (mut info, headers) = split_and_match(&header_re, &ansi_string);

    // this should find any ANSI formatted clap argument in the info lines
    let arg_re = Regex::new(&format!(r"(?m)^\s*{BOLD_RE}.*{RESET_RE}")).unwrap();

    for info_lines in info
        .iter_mut()
        // skip three sets of info lines,
        // since the "about" section
        // is before the first header
        .skip(3)
    {
        let (helps, args) = split_and_match(&arg_re, info_lines);

        *info_lines = args
            .into_iter()
            .zip(helps.into_iter().skip(1))
            .enumerate()
            .map(|(i, (arg, help))| {
                let color = if i % 2 == 0 { BRIGHT_WHITE } else { WHITE };
                format!("{color}{arg}{color}{help}{RESET}")
            })
            .collect();
    }

    // this prints "about"
    print!("{}", info[0]);

    for (h, i) in headers.into_iter().zip(info.into_iter().skip(1)) {
        print!("{h}{i}");
    }
}
