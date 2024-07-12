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

### Changed

### Deprecated

### Removed

### Fixed


## [0.2.0] - 2024-07-12

### Added

- added Nats and Bits structs
- added AlignmentBuilder struct
- added Profile::raise_relative_entropy()
- added Profile::calibrate_tau()
- added Hmm::from_blosum_62_and_sequence()
- added Sequence::random_amino()
- added several AntiDiagonal methods:
    - valid()
    - reset()
    - is_default()
    - is_single_cell()
    - intersects()
    - merge() 
    - replace_with()
    - grow_up()
    - grow_down()
    - grow_left()
    - grow_right()
- added several AntiDiagonalBounds methods:
    - square_corners()
    - advance_forward()
    - advance_reverse()
    - bounding_box()
    - left_offset_to()
    - right_offset_to()
    - cloud_relationship()
    - vec_image()
    - ascii()
- added many unit tests for AntiDiagonal and AntiDiagonalBounds
- added unit tests for DpMatrixFlat and DpMatrixSparse
- added `debug` Cargo feature:
    - image crate dependency
    - AntiDiagonalBounds::image() 
- added VecUtils trait for some common numerical vector operations

### Changed

- Nats and Bits structs are now used to represent alignment scores
- most fields on the Alignment struct are now Option<T>
- renamed CloudBound to AntiDiagonal
- renamed AntiDiagonalBounds::join_merge() to merge() and improved merging logic to handle more edge cases
- added score field to Seed struct
- refactored and improved tabular output

### Removed

- removed DPMatrix3D struct

