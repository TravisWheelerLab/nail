mod align;
use std::{
    collections::HashMap,
    io::{stdout, Write},
    sync::{Arc, Mutex},
};

pub use align::*;
mod prep;
use libnail::{
    align::{
        backward, cloud_search_backward, cloud_search_forward, forward, null_one_score,
        null_two_score, optimal_accuracy, posterior,
        structs::{
            Alignment, AlignmentBuilder, AntiDiagonalBounds, CloudMatrixLinear, DpMatrixSparse,
            RowBounds, Seed, Trace,
        },
        traceback, CloudSearchParams, CloudSearchScores,
    },
    output::output_tabular::TableFormat,
    structs::{Profile, Sequence},
};
pub use prep::*;
mod search;
use rayon::iter::{IntoParallelIterator, ParallelIterator};
pub use search::*;
mod seed;
pub use seed::*;

use crate::{
    database::{Database, ProfileCollection, SequenceCollection},
    extension_traits::PathBufExt,
};

#[derive(Default, Clone)]
struct SeedStep {
    seeds: HashMap<String, HashMap<String, Seed>>,
}

impl SeedStep {
    fn run(&self, profile: &Profile, target: &Sequence) -> Option<Seed> {
        Some(self.seeds.get(&profile.name)?.get(&target.name)?.clone())
    }
}

#[derive(Default, Clone)]
struct CloudSearchStep {
    cloud_matrix: CloudMatrixLinear,
    forward_bounds: AntiDiagonalBounds,
    reverse_bounds: AntiDiagonalBounds,
    row_bounds: RowBounds,
    params: CloudSearchParams,
}

impl CloudSearchStep {
    // TODO: move this to libnail
    fn cloud_score(forward_scores: &CloudSearchScores, reverse_scores: &CloudSearchScores) -> f32 {
        // this approximates the score for the forward
        // cloud that extends past the seed end point
        let disjoint_forward_score = forward_scores.max_score - forward_scores.max_score_within;

        // this approximates the score for the reverse
        // cloud that extends past the seed start point
        let disjoint_reverse_score = reverse_scores.max_score - reverse_scores.max_score_within;

        // this approximates the score of the intersection
        // of the forward and reverse clouds
        let intersection_score = forward_scores
            .max_score_within
            .max(reverse_scores.max_score_within);

        let cloud_score_nats = intersection_score + disjoint_forward_score + disjoint_reverse_score;

        // convert to bits
        cloud_score_nats / std::f32::consts::LN_2
    }

    fn run(&mut self, profile: &Profile, target: &Sequence, seed: &Seed) -> Option<&RowBounds> {
        self.cloud_matrix.reuse(profile.length);
        self.forward_bounds.reuse(target.length, profile.length);
        self.reverse_bounds.reuse(target.length, profile.length);
        self.row_bounds.reuse(target.length);

        let forward_scores = cloud_search_forward(
            profile,
            target,
            seed,
            &mut self.cloud_matrix,
            &self.params,
            &mut self.forward_bounds,
        );

        let reverse_scores = cloud_search_backward(
            profile,
            target,
            seed,
            &mut self.cloud_matrix,
            &self.params,
            &mut self.reverse_bounds,
        );

        // TODO: make this actually work
        // if Self::cloud_score(&forward_scores, &reverse_scores) < self.params.threshold {
        //     return None;
        // }
        //
        // let cloud_pvalue = (-profile.forward_lambda as f64
        //     * (cloud_score_bits as f64 - profile.forward_tau as f64))
        //     .exp();

        // if cloud_pvalue >= data.cloud_pvalue_threshold {
        //     continue;
        // }

        let bounds_intersect =
            self.forward_bounds.max_anti_diagonal_idx >= self.reverse_bounds.min_anti_diagonal_idx;

        if bounds_intersect {
            AntiDiagonalBounds::join_merge(&mut self.forward_bounds, &self.reverse_bounds);
            if !self.forward_bounds.valid() {
                self.forward_bounds.fill_rectangle(
                    seed.target_start,
                    seed.profile_start,
                    seed.target_end,
                    seed.profile_end,
                );
            }

            self.forward_bounds.square_corners();
            self.forward_bounds.trim_wings();

            self.row_bounds
                .fill_from_anti_diagonal_bounds(&self.forward_bounds);

            if !self.row_bounds.valid() {
                self.row_bounds.fill_rectangle(
                    seed.target_start,
                    seed.profile_start,
                    seed.target_end,
                    seed.profile_end,
                );
            }
        } else {
            self.row_bounds.fill_rectangle(
                seed.target_start,
                seed.profile_start,
                seed.target_end,
                seed.profile_end,
            );
        }

        Some(&self.row_bounds)
    }
}

#[derive(Default, Clone)]
struct AlignmentStep {
    forward_matrix: DpMatrixSparse,
    backward_matrix: DpMatrixSparse,
    posterior_matrix: DpMatrixSparse,
    optimal_matrix: DpMatrixSparse,
    forward_pvalue_threshold: f64,
    evalue_threshold: f64,
}

impl AlignmentStep {
    fn run(
        &mut self,
        profile: &mut Profile,
        target: &Sequence,
        bounds: &RowBounds,
    ) -> Option<Alignment> {
        // configuring for the target length
        // adjusts special state transitions
        profile.configure_for_target_length(target.length);

        self.forward_matrix
            .reuse(target.length, profile.length, bounds);
        self.backward_matrix
            .reuse(target.length, profile.length, bounds);
        self.posterior_matrix
            .reuse(target.length, profile.length, bounds);
        self.optimal_matrix
            .reuse(target.length, profile.length, bounds);

        // we use the forward score to compute the final bit score (later)
        let forward_score_nats = forward(profile, target, &mut self.forward_matrix, bounds);

        // for now we compute the P-value for filtering purposes
        let forward_pvalue = (-profile.forward_lambda as f64
            * ((forward_score_nats / std::f32::consts::LN_2) as f64 - profile.forward_tau as f64))
            .exp();

        // if forward_pvalue >= self.forward_pvalue_threshold {
        //     return None;
        // }

        backward(profile, target, &mut self.backward_matrix, bounds);

        posterior(
            profile,
            &self.forward_matrix,
            &self.backward_matrix,
            &mut self.posterior_matrix,
            bounds,
        );

        optimal_accuracy(
            profile,
            &self.posterior_matrix,
            &mut self.optimal_matrix,
            bounds,
        );

        let mut trace = Trace::new(target.length, profile.length);
        traceback(
            profile,
            &self.posterior_matrix,
            &self.optimal_matrix,
            &mut trace,
            bounds.target_end,
        );

        let null_one = null_one_score(target.length);
        let null_two = null_two_score(&self.posterior_matrix, profile, target, bounds);

        match AlignmentBuilder::new(&trace)
            .with_profile(profile)
            .with_target(target)
            .with_target_count(1)
            .with_raw_score(forward_score_nats)
            .with_null_one(null_one)
            .with_null_two(null_two)
            .with_cell_fraction(bounds.num_cells() as f32 / (profile.length * target.length) as f32)
            .build()
        {
            Ok(alignment) => Some(alignment),
            Err(_) => None,
        }
    }
}

#[derive(Default, Clone)]
pub struct Pipeline {
    seed: SeedStep,
    cloud_search: CloudSearchStep,
    align: AlignmentStep,
}

impl Pipeline {
    fn run(&mut self, profile: &mut Profile, target: &Sequence) -> Option<Alignment> {
        self.seed
            .run(profile, target)
            .and_then(|seed| self.cloud_search.run(profile, target, &seed))
            .and_then(|bounds| self.align.run(profile, target, bounds))
    }
}

pub fn search_new(args: &SearchArgs) -> anyhow::Result<()> {
    {
        // quickly make sure we can write the results
        args.output_args.tsv_results_path.open(true)?;
    }

    let seeds_path = args.prep_dir.path.join("./seeds.json");
    let prep_args = PrepArgs {
        query_path: args.query_path.clone(),
        target_path: args.target_path.clone(),
        skip_hmmbuild: args.prebuilt_query_hmm_path.is_some(),
        prep_dir: args.prep_dir.clone(),
        common_args: args.common_args.clone(),
    };

    let seed_args = SeedArgs {
        prep_dir_path: Default::default(),
        seeds_path: seeds_path.clone(),
        prebuilt_query_hmm_path: args.prebuilt_query_hmm_path.clone(),
        prep_dir: args.prep_dir.clone(),
        mmseqs_args: args.mmseqs_args.clone(),
        common_args: args.common_args.clone(),
    };

    let align_args = AlignArgs {
        query_path: args.query_path.clone(),
        target_path: args.target_path.clone(),
        seeds_path,
        nail_args: args.nail_args.clone(),
        output_args: args.output_args.clone(),
        common_args: args.common_args.clone(),
    };

    prep(&prep_args)?;
    let (profiles_vec, seed_map) = seed(&seed_args)?;

    let mut seeds: HashMap<String, HashMap<String, Seed>> = HashMap::new();
    seed_map.iter().for_each(|(profile, seed_vec)| {
        let profile_map = seeds.entry(profile.to_string()).or_default();
        seed_vec.iter().for_each(|seed| {
            profile_map.insert(seed.target_name.clone(), seed.clone());
        });
    });

    let profiles = ProfileCollection::new(profiles_vec);
    let seqs = SequenceCollection::new(Sequence::amino_from_fasta(align_args.target_path)?);

    let pipeline = Pipeline {
        seed: SeedStep { seeds },
        ..Default::default()
    };

    let output = Arc::new(Mutex::new(Output {
        alignment_writer: Box::new(stdout()),
        table_writer: Box::new(args.output_args.tsv_results_path.open(true)?),
        table_format: TableFormat::new(&DEFAULT_COLUMNS)?,
        header_status: HeaderStatus::Unwritten,
    }));

    run_pipeline(&profiles, &seqs, pipeline, output);

    Ok(())
}

pub enum HeaderStatus {
    Unwritten,
    Written,
}

pub struct Output {
    alignment_writer: Box<dyn Write + Send>,
    table_writer: Box<dyn Write + Send>,
    table_format: TableFormat,
    header_status: HeaderStatus,
}

impl Output {
    pub fn write(&mut self, alignments: &mut [Alignment]) -> anyhow::Result<()> {
        self.table_format.update_widths(alignments);

        alignments.sort_by(|a, b| a.e_value.partial_cmp(&b.e_value).unwrap());

        if let HeaderStatus::Unwritten = self.header_status {
            let header = TableFormat::header(&self.table_format)?;
            writeln!(self.table_writer, "{header}")?;
            self.header_status = HeaderStatus::Written;
        }

        alignments.iter().for_each(|ali| {
            writeln!(
                self.table_writer,
                "{}",
                ali.tab_string_formatted(&self.table_format)
            )
            .expect("failed to write tabular output");

            writeln!(self.alignment_writer, "{}", ali.ali_string())
                .expect("failed to write alignment output");
        });

        Ok(())
    }
}

pub fn run_pipeline<A, B>(
    profiles: &dyn Database<A, Arc<Mutex<Profile>>>,
    targets: &(dyn Database<B, Arc<Sequence>> + Send + Sync),
    pipeline: Pipeline,
    output: Arc<Mutex<Output>>,
) where
    A: Sync + Send,
    B: Sync + Send,
{
    profiles.into_par_iter().for_each_with(
        (pipeline, targets, output),
        |(pipeline, targets, output), profile| {
            let mut profile_guard = profile.lock().unwrap();
            let mut alignments: Vec<Alignment> = targets
                .iter()
                .filter_map(|target| pipeline.run(&mut profile_guard, &target))
                .collect();

            match output.lock() {
                Ok(mut guard) => {
                    guard.write(&mut alignments).unwrap();
                }
                Err(_) => panic!(),
            }
        },
    )
}
