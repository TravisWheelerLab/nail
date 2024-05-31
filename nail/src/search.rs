use std::collections::HashMap;
use std::process::Command;
use std::sync::{Arc, Mutex};

use crate::args::{AlignArgs, SearchArgs};
use crate::pipeline::{
    run_pipeline_profile_to_sequence, run_pipeline_sequence_to_sequence, DefaultAlignStep,
    DefaultCloudSearchStep, DefaultSeedStep, DoubleSeedStep, FullDpCloudSearchStep, OutputStep,
    Pipeline,
};
use crate::util::{
    check_hmmer_installed, guess_query_format_from_query_file, CommandExt, FileFormat, PathBufExt,
};

use libnail::structs::hmm::parse_hmms_from_p7hmm_file;
use libnail::structs::{Profile, Sequence};

use anyhow::Context;

pub enum Queries {
    Sequence(Vec<Sequence>),
    Profile(Vec<Profile>),
    DoubleProfile((Vec<Profile>, Vec<Profile>)),
}

pub fn search(args: &SearchArgs) -> anyhow::Result<()> {
    {
        // quickly make sure we can write the results
        args.output_args.tsv_results_path.open(true)?;
        if let Some(path) = &args.output_args.ali_results_path {
            path.open(true)?;
        }
    }

    let seeds_path = args.prep_dir.join("./seeds.json");

    let align_args = AlignArgs {
        query_path: args.query_path.clone(),
        target_path: args.target_path.clone(),
        seeds_path,
        nail_args: args.nail_args.clone(),
        output_args: args.output_args.clone(),
        common_args: args.common_args.clone(),
    };

    let query_format = guess_query_format_from_query_file(&args.query_path)?;

    let queries = match query_format {
        FileFormat::Fasta => {
            let queries = Sequence::amino_from_fasta(&args.query_path)
                .context("failed to read query fasta")?;

            Queries::Sequence(queries)
        }
        FileFormat::Hmm => {
            let queries: Vec<Profile> = parse_hmms_from_p7hmm_file(&args.query_path)
                .context("failed to read query hmm")?
                .iter()
                .map(Profile::new)
                .collect();
            Queries::Profile(queries)
        }
        FileFormat::Stockholm => {
            check_hmmer_installed()?;
            let hmm_path = args.prep_dir.join("query.hmm");

            Command::new("hmmbuild")
                .args(["--cpu", &args.common_args.num_threads.to_string()])
                .arg(&hmm_path)
                .arg(&args.query_path)
                .run()?;

            let queries_a = parse_hmms_from_p7hmm_file(&hmm_path)?
                .iter()
                .map(Profile::new)
                .collect();

            Command::new("hmmbuild")
                .args(["--ere", "1.0"])
                .args(["--cpu", &args.common_args.num_threads.to_string()])
                .arg(&hmm_path)
                .arg(&args.query_path)
                .run()?;

            let queries_b = parse_hmms_from_p7hmm_file(&hmm_path)?
                .iter()
                .map(Profile::new)
                .collect();

            Queries::DoubleProfile((queries_a, queries_b))
        }
        _ => {
            panic!()
        }
    };

    let targets = Sequence::amino_from_fasta(&align_args.target_path)
        .context("failed to read target fasta")?;

    let pipeline = Pipeline {
        seed: match queries {
            Queries::DoubleProfile(_) => Box::new(DoubleSeedStep::new(
                &queries,
                &targets,
                &args.prep_dir,
                args.common_args.num_threads,
                &args.mmseqs_args,
            )?),
            _ => Box::new(DefaultSeedStep::new(
                &queries,
                &targets,
                &args.prep_dir,
                args.common_args.num_threads,
                &args.mmseqs_args,
            )?),
        },

        cloud_search: match args.nail_args.full_dp {
            true => Box::<FullDpCloudSearchStep>::default(),
            false => Box::new(DefaultCloudSearchStep::new(&align_args)),
        },
        align: Box::new(DefaultAlignStep::new(&align_args, targets.len())),
        output: Arc::new(Mutex::new(OutputStep::new(&args.output_args)?)),
    };

    let target_map: HashMap<String, Sequence> =
        targets.into_iter().map(|t| (t.name.clone(), t)).collect();

    match queries {
        Queries::Sequence(queries) => {
            run_pipeline_sequence_to_sequence(&queries, &target_map, pipeline);
        }
        Queries::Profile(mut queries) => {
            run_pipeline_profile_to_sequence(&mut queries, &target_map, pipeline);
        }
        Queries::DoubleProfile((mut queries, _)) => {
            run_pipeline_profile_to_sequence(&mut queries, &target_map, pipeline);
        }
    }

    Ok(())
}
