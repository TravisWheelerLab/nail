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

### Added
- added support for `--double-seed` for sequence to sequence search
- added `Fasta.names_iter()`
- added `mmseqs::consts::{BLOSUM_62, BLOSUM_80}`
- added `MmseqsScoreModel`
- added `MmseqsDbPaths.dir()`

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

