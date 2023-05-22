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
- Alignment struct now stores P-values
- Added null score adjustment
- Added CloudBoundGroup.json(), CloudBoundGroup.valid() and CloudBoundGroup.bounds()
- Added RowBounds.json()
- Added write_standard_output() 

### Changed
- Tabular output header lines now start with "#"
- Replaced UnknownSequenceCharacterError with UnknownUtf8SequenceByteError and UnknownDigitalSequenceByteError
- Alignment::new() now takes a ScoreParams struct
- forward_bounded() now returns the score in the last C state (target_end)
- Profile, CloudMatrixLinear, DpMatrixSparse now derive Clone
- E-Values are now f64 values instead of f32
- Renamed P7Alphabet to Alphabet

### Removed
- Removed CloudDebugAnnotations struct

### Fixed
- Fixed E-Value column width in write_tabular_output()
- Fixed parameter names on DpMatrix.get_special() and DpMatrix.set_special()
- Fixed float parsing regex to capture negative signs
- Fixed fasta header parsing 
- Fixed a traceback bug that caused erroneous J state usage

# [0.1.1] - 2023-04-16

### Changed
- Added E-value calculation to the Alignment struct
- Added "bit score" field to tabular output

## [0.1.0] - 2023-04-15

- Initial release

