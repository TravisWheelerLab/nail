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


## [0.4.0] - 2025-6-18

### Added
- added `PartialEq` and `PartialOrd` impls for `Bits` and `Nats`
- added struct `Cell`
- added struct `CloudImage` for better cloud visualization API
- added structs `AdMatrixSparse`, `AdMatrixLinear`, `AdMatrixQuadratic`
- added enum `AminoAcid`
- added traits `CoreCellIndexable`, `BackgroundCellIndexable`
- added enums `CoreState`, `BackgroundState`
- added structs `Emission`, `CoreToCore`, `BackgroundLoop`, `CoreEntry`
- added traits `MinAssign`, `MaxAssign` and blanket impls
- added macro `assert_eq_pairs!`
### Changed
- `Profile.match_scores` and `Profile.insert_scores` vectors combined into `Profile.emission_scores`
- structs `Sequence`, `Profile` now have an ending pad position
- changed field names on `Seed` struct to use "seq" and "prf" instead of "target" and "profile"
- renamed struct `AntiDiagonal` to `Bound`
- renamed struct `AntiDiagonalBounds` to `Cloud`
- changed cloud search implementation to better match the Forward/Backward recurrence
- renamed trait `VecUtils` to `VecMath`
### Removed
- removed structs `CloudAntiDiagonal`, `CloudMatrixLinear`

## [0.3.0] - 2024-12-23

### Added
- added `min()` to `Score`, `Nats`, and `Bits`
- added `Alignment::vert_string()`
- added consts `Alignment::PROFILE_GAP_BYTE` and `Alignment::TARGET_GAP_BYTE`
- added `TraceStep::is_core_state()`

### Changed
- `Sequence` and `DpMatrixSparse` and now call `shrink_to_fit()` on internal data
- `Profile::raise_relative_entropy()` is now `adjust_mean_relative_entropy()`
    - the mean relative entropy (MRE) can now be raised or lowered
    - the algorithm for raising (MRE) has been improved
- `util::avg_relative_entropy()` is now `mean_relative_entropy()`
- parse_hmms_from_p7hmm_file(path) is now Hmm::from_p7hmm(buf)
- refactored Alignment struct
    - fields are now grouped into sub-structs
    - refactored AlignmentBuilder for changes
    - added impl for `AsRef\<Alignment\>`
- `DpMatrixSparse` should now resize properly
- `CloudSearchScores` renamed to `CloudSearchResults`, added `num_cells_computed` field
- `TableFormat::update_widths()` now takes a slice of `AsRef\<Alignment\>` instead of `&\[Alignment\]`
- changed the formatting of the output of `Alignment::ali_string()`

### Fixed
- Fixed the score returned from `Forward()` so that it includes N/C state
  transitions from outside of the cloud
- Fixed an alignment boundary computation bug


## [0.2.0] - 2024-07-12

### Added

- added `Nats` and `Bits` structs
- added `AlignmentBuilder` struct
- added `Profile::raise_relative_entropy()`
- added `Profile::calibrate_tau()`
- added `Hmm::from_blosum_62_and_sequence()`
- added `Sequence::random_amino()`
- added several `AntiDiagonal` methods:
    - `valid()`
    - `reset()`
    - `is_default()`
    - `is_single_cell()`
    - `intersects()`
    - `merge() `
    - `replace_with()`
    - `grow_up()`
    - `grow_down()`
    - `grow_left()`
    - `grow_right()`
- added several `AntiDiagonalBounds` methods:
    - `square_corners()`
    - `advance_forward()`
    - `advance_reverse()`
    - `bounding_box()`
    - `left_offset_to()`
    - `right_offset_to()`
    - `cloud_relationship()`
    - `vec_image()`
    - `ascii()`
- added many unit tests for `AntiDiagonal` and `AntiDiagonalBounds`
- added unit tests for `DpMatrixFlat` and `DpMatrixSparse`
- added `debug` Cargo feature:
    - `image` crate dependency
    - `AntiDiagonalBounds::image()` 
- added `VecUtils` trait for some common numerical vector operations

### Changed

- `Nats` and `Bits` structs are now used to represent alignment scores
- most fields on the `Alignment` struct are now `Option<T>`
- renamed `CloudBound` to `AntiDiagonal`
- renamed `AntiDiagonalBounds::join_merge()` to `merge()` and improved merging logic to handle more edge cases
- added `score` field to `Seed` struct
- refactored and improved tabular output

### Removed

- removed `DPMatrix3D` struct

