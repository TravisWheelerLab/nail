# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

<!---
### Added
### Changed
### Deprecated
### Removed
### Fixed
### Security
-->


## [Unreleased]

### Added

- added Stats module for keeping track of pipeline summary statistics

### Changed

### Deprecated

### Removed

### Fixed


## [0.2.0] - 2024-07-12

### Added

- added support for sequence to sequence search
- added support for HMM query files
- added Pipeline struct
- added SeedStep, CloudStep, AlignStep, and OutputStep traits
- added MmseqsArgs struct
- added write_mmseqs_database() function
- added write_mmseqs_profile_database() function
- added run_mmseqs_search() function
- added seeds_from_mmseqs_align_tsv() function


### Changed

- large scale refactor of file structure
    - all pipeline code was moved around to various files
    - args structs better consolidated in args.rs
    - extension traits moved to util.rs
- no longer depend on `hmmbuild` at runtime
- the alignment pipeline is now flexibly configured via the Pipeline struct and step traits
- the default pipeline now runs `mmseqs search` twice at two profile relative entropy targets

### Removed

- removed `seed` and `prep` commands from the CLI
- removed check_hmmer_installed()
- removed viz.rs

