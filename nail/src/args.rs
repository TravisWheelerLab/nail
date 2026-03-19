use std::{io::Write, path::PathBuf, str::FromStr};

use anyhow::bail;
use clap::{Args, Parser, Subcommand};

use crate::util::{term::*, PathExt};

#[derive(Subcommand)]
#[allow(clippy::large_enum_variant)]
pub enum NailSubCommands {
    #[command(about = "Run nail's protein search pipeline")]
    Search(SearchArgs),
    #[command(subcommand, hide = true)]
    Dev(DevSubCommands),
}

#[derive(Subcommand)]
#[allow(clippy::large_enum_variant)]
pub enum DevSubCommands {
    #[command()]
    Play(SearchArgs),
    Search(SearchArgs),
    Mx(SearchArgs),
}

#[derive(Parser)]
#[command(version)]
#[command(name = "nail")]
#[command(
    about = "Using MMseqs2 to find rough alignment seeds, perform bounded profile HMM sequence alignment"
)]
pub struct NailCli {
    #[command(subcommand)]
    pub command: NailSubCommands,
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

    #[command(flatten)]
    #[clap(next_help_heading = "File I/O options")]
    pub io_args: IoArgs,

    #[command(flatten)]
    #[clap(next_help_heading = "Pipeline options")]
    pub pipeline_args: PipelineArgs,

    /// Arguments that are passed to MMseqs2
    #[command(flatten)]
    #[clap(next_help_heading = "Seeding options")]
    pub mmseqs_args: MmseqsArgs,

    #[command(flatten)]
    #[clap(next_help_heading = "Expert options")]
    pub expert_args: ExpertArgs,

    #[command(flatten)]
    #[clap(next_help_heading = "Dev options")]
    pub dev_args: DevArgs,
}

impl SearchArgs {
    pub fn validate(&mut self) -> anyhow::Result<()> {
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

        {
            // quickly make sure we can write to all of the results paths
            self.io_args.temp_dir_path.create_dir()?;

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

        if self.mmseqs_args.prog_seed {
            if self.mmseqs_args.max_seqs != 2_147_483_647 {
                bail!(
                    "the argument '{YELLOW}--mmseqs-max-seqs{RESET}' is set wrong: {}",
                    self.mmseqs_args.max_seqs
                )
            }

            if self.mmseqs_args.prog_n.is_none() {
                bail!("the argument '{YELLOW}--prog-n{RESET}' is unset")
            }

            if self.mmseqs_args.prog_n.is_none() {
                bail!("the argument '{YELLOW}--prog-f{RESET}' is unset")
            }
        } else {
            #[allow(clippy::collapsible_if)]
            if self.mmseqs_args.prog_n.is_some() {
                bail!("the argument '{YELLOW}--prog-n{RESET}' cannot be used without '{YELLOW}--prog-seed{RESET}'")
            }
            if self.mmseqs_args.prog_n.is_some() {
                bail!("the argument '{YELLOW}--prog-f{RESET}' cannot be used without '{YELLOW}--prog-seed{RESET}'")
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
        writeln!(out, " ├─ mmseqs -k: {}", self.mmseqs_args.k)?;
        writeln!(out, " ├─ mmseqs -s: {:.}", self.mmseqs_args.s)?;
        writeln!(
            out,
            " ├─ mmseqs --max-seqs: {:.}",
            self.mmseqs_args.max_seqs
        )?;
        writeln!(out, " ├─ prog-seed: {}", self.mmseqs_args.prog_seed)?;
        if self.mmseqs_args.prog_seed {
            writeln!(
                out,
                "   ├─ n: {}",
                self.mmseqs_args.prog_n.unwrap_or_default()
            )?;
            writeln!(
                out,
                "   └─ f: {}",
                self.mmseqs_args.prog_f.unwrap_or_default()
            )?;
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
    #[arg(long = "tmp-dir", default_value = "tmp-nail/", value_name = "PATH")]
    pub temp_dir_path: PathBuf,

    /// Allow nail to overwrite files
    #[arg(long = "allow-overwrite", default_value_t = false)]
    pub allow_overwrite: bool,
}

#[derive(Args, Debug, Clone, Default)]
pub struct PipelineArgs {
    /// Pruning parameter alpha
    #[arg(
        short = 'A',
        default_value_t = 10.0,
        value_name = "X",
        help = "Cloud search parameter α:\n  \
                local score pruning threshold"
    )]
    pub alpha: f32,

    /// Pruning parameter beta
    #[arg(
        short = 'B',
        default_value_t = 16.0,
        value_name = "X",
        help = "Cloud search parameter β:\n  \
                global score pruning threshold"
    )]
    pub beta: f32,

    /// Pruning parameter gamma
    #[arg(
        short = 'G',
        default_value_t = 5,
        value_name = "N",
        help = "Cloud search parameter γ:\n  \
                at minimum, compute N anti-diagonals"
    )]
    pub gamma: usize,

    /// Seeding filter threshold
    #[arg(
        short = 'S',
        default_value_t = 1e-4,
        value_name = "X",
        help = "Seeding filter threshold:\n  \
                filter hits with P-value > X"
    )]
    pub seed_pvalue_threshold: f64,

    /// Cloud search filter threshold
    #[arg(
        short = 'C',
        default_value_t = 1e-2,
        value_name = "X",
        help = "Cloud search threshold:\n  \
                filter hits with P-value > X"
    )]
    pub cloud_pvalue_threshold: f64,

    /// Forward filter threshold
    #[arg(
        short = 'F',
        default_value_t = 1e-4,
        value_name = "X",
        help = "Forward filter threshold:\n  \
                filter hits with P-value > X"
    )]
    pub forward_pvalue_threshold: f64,

    /// Final E-value threshold
    #[arg(
        short = 'E',
        default_value_t = 10.0,
        value_name = "X",
        help = "Final reporting threshold:\n  \
                filter hits with E-value > X"
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
pub struct MmseqsArgs {
    /// MMseqs2 Parameter: k-mer length (0: automatically set to optimum)
    #[arg(
        long = "mmseqs-k",
        default_value_t = 6usize,
        value_name = "N",
        help = "MMseqs2 parameter:\n  \
                k-mer length (0: automatically set to optimum)"
    )]
    pub k: usize,

    /// MMseqs2 Parameter: Sensitivity: 1.0 faster; 4.0 fast; 7.5 sensitive
    #[arg(
        long = "mmseqs-s",
        default_value_t = 10.0,
        value_name = "X",
        help = "MMseqs2 parameter:\n  \
                Sensitivity: 1.0 faster; 4.0 fast; 7.5 sensitive"
    )]
    pub s: f32,

    /// MMseqs2 Parameter: Maximum results per query sequence allowed to pass the prefilter
    #[arg(
        long = "mmseqs-max-seqs",
        default_value_t = 200usize,
        default_value_if("prog_seed", "true", "2147483647"),
        value_name = "N",
        help = "MMseqs2 parameter:\n  \
                Maximum results per query sequence allowed to pass the prefilter"
    )]
    pub max_seqs: usize,

    /// MMseqs2 Parameter: Correct for locally biased amino acid composition (range 0-1)
    #[arg(
        long = "mmseqs-comp-bias-corr",
        default_value = None,
        value_name = "N",
        hide = true,
        help = "MMseqs2 parameter:\n  \
                Correct for locally biased amino acid composition (range 0-1)"
    )]
    pub comp_bias_corr: Option<usize>,

    /// Enable progressive seeding
    #[arg(long, action, conflicts_with = "max_seqs")]
    pub prog_seed: bool,

    /// The initial number of mmseqs alignments per query
    #[arg(
        long,
        default_value = "200",
        default_value_if("prog_seed", "true", "200"),
        default_value_if("prog_seed", "false", None),
        value_name = "N",
        conflicts_with = "max_seqs",
        help = "Progressive seeding:\n  \
                the initial number of mmseqs alignments per query"
    )]
    pub prog_n: Option<usize>,

    /// test
    #[arg(
        long,
        default_value = "0.01",
        default_value_if("prog_seed", "true", "0.01"),
        default_value_if("prog_seed", "false", None),
        value_name = "X",
        conflicts_with = "max_seqs",
        help = "Progressive seeding:\n  \
                the fraction of hits required to continue progressive seeding"
    )]
    pub prog_f: Option<f32>,
}
