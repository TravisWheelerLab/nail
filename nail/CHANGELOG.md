# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

<!---
## [Unreleased]
### Added
### Changed
### Deprecated
### Removed
### Fixed
### Security
-->

## [Unreleased]


## [0.5.0] - 2026-3-19

### Added
- added `nail dev` command to CLI
- added `Fasta::from_path_par`
- added `LexicalFastaIndex::from_path`
- added `FileExt` trait for `File` API extensions
- added `Fasta.lengths_iter()`
- added `write_seed_map()`
- added `read_seed_map()`
- added mod `io::database`
- added mod `io::index`
- added mod `io::p7hmm`
- added mod `io::seeds`
- added `SearchArgs::{validate(), write()}`
- added CLI params: `--prog-seed`, `--prog-n`, `--prog-f`
- added consts `mmseqs::consts::{ALIGN_DBTYPE, PREFILTER_DBTYPE}`
- added structs `mmseqs::PrefilterDb, Descriptor, ByteBuffer`
- added methods `MmseqsDbPaths:{destroy(), check(), db_exists(), rmeove_db()}`
- added function `mmseqs:{run_mmseqs_convertalis()}`
- added trait `TableDisplay`
- added method `PipelineResult::stat_string()`
- added enum `stats::MmseqsTimed`
- added mod `util::term`
- added dev CLi options to adjust formatting numeric precision

### Changed
- changed `FastaOffset` length fields to include `_bytes` suffix
- mod `io::rayon` renamed to `io::impl_rayon`
- `LexicalFastaIndex::new()` now also takes `start: Option<u64>` to supply a relative starting offset
- moved `Fasta` struct & related to mod `io::fasta`
- changed mmseqs2 CLI params:
    - removed `--mmseqs-k-score` and `--min-ungapped-score`
    - added `--mmseqs-s`
    - changed default `--mmseqs-k` parameter to k=6
- changed `IoArgs::temp_dir_path` default from "/tmp" to "/tmp-nail"
- split function `run_mmseqs_search()` into `run_mmseqs_align()` and `run_mmseqs_prefilter()`
- added functions `seed_max_seqs(), seed_progressive()`
- changed search pipeline:
    - SeedStage is no longer used (for now)
    - added fields `Pipeline:{profiles, prf}`
    - pipeline now parallelizes over seeds instead of profiles (better performance for unbalanced seed distributions)

### Removed
- removed `--double-seed` option from CLI
- removed method `tab_string()` from structs `AlignStageResult`, CloudStageResult`, `PipelineResult`
- removed unused error structs `ProfileNotFoundError`, `TargetNotFoundError`
- removed field `Pipeline:seed`
- removed struct `DefaultSeedStage`
- removed functions `seed_profile_to_sequence(), seed_sequence_to_sequence()`

### Fixed
- fixed a bug in `FastaParser::offset()`

## [0.4.0] - 2025-6-18

### Added
- added `--version` commmand to CLI
- added `--max-seed` command to dev CLI
- added `MaxSeedStage` to force a seed for every alignment
- added support for `--double-seed` for sequence to sequence search
- added `Fasta.names_iter()`
- added `mmseqs::consts::{BLOSUM_62, BLOSUM_80}`
- added `MmseqsScoreModel`
- added `MmseqsDbPaths.dir()`

### Changed
- changed default `alpha` and `beta` parameters to `10.0` and `16.0`

## [0.3.0] - 2024-12-23

### Added

- added module `io`
    - `Fasta` struct for `Sequence` retrieval streamed from disk
    - `SequenceDatabase` trait for generic search databases
    - submodule `io::rayon`, which includes impls for `rayon` parallel iterators over `SequenceDatabase`s
- added jemalloc support with feature flag `jemalloc`
- added `Stats` module for keeping track of pipeline summary statistics
- search pipeline now uses the `thread_local` crate instead of letting Rayon clone whenever it wants to
- added `stats_results_path` to `OutputArgs` struct
- added `StageResult<S, D>` enum where:
    - `StageResult::Filtered` indicates that the result was filtered, and
    - `StageResult::Passed` indicates that the result should be passed to the next stage
- added `PipelineResult` struct, which captures the `StageResult` returned by each pipeline stage
- added dependency on the `derive_builder` crate
- added `CloudStageStats`, `AlignStageStats`, and `OutputStageStats` structs, each of which has a derived builder struct from the `derive_builder` crate
- added dependency on `indexmap` crate


### Changed

- `OutputArgs.evalue_threshold` field renamed to `e_value_threshold`
- refactored `Pipeline` struct:
    - `output` field is now a clonable `OutputStep` instead of `Arc<Mutex<OutputStep>>`
    - `run()` now returns `anyhow::Result<()>` instead of `anyhow::Result<Vec<Alignment>>`
    - `run()` now keeps better track of filtered results and relevant stats using the `PipelineResult` struct
- refactored pipeline stage structs:
    - they have been renamed from Step to Stage
    - each stage now returns a `StageResult<S, D>` where:
        - `S` is a struct containing observations for pipeline summary statistics
        - `D` is a struct containing the data structure relevant to the stage
    - `OutputStage` now:
        - holds an `Arc<Mutex<T>>` for each individual writer
        - can optionally write any or none of the output categories
        - can now write summary statistics to the `stats_results_path` data for each call to `Pipeline.run()`
- refactored CLI


## [0.2.0] - 2024-07-12

### Added

- added support for sequence to sequence search
- added support for HMM query files
- added `Pipeline` struct
- added `SeedStep`, `CloudStep`, `AlignStep`, and `OutputStep` traits
- added `MmseqsArgs` struct
- added `write_mmseqs_database()` function
- added `write_mmseqs_profile_database()` function
- added `run_mmseqs_search()` function
- added `seeds_from_mmseqs_align_tsv()` function


### Changed

- large scale refactor of file structure
    - all pipeline code was moved around to various files
    - args structs better consolidated in args.rs
    - extension traits moved to `util.rs`
- no longer depend on `hmmbuild` at runtime
- the alignment pipeline is now flexibly configured via the Pipeline struct and step traits
- the default pipeline now runs `mmseqs search` twice at two profile relative entropy targets

### Removed

- removed `seed` and `prep` commands from the CLI
- removed `check_hmmer_installed()`
- removed `viz.rs`

