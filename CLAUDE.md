# nail

A profile HMM sequence aligner, amino acids only. `nail` is the command-line
tool; `libnail` is the library holding the sparse alignment algorithms. nail
runs MMseqs2 `search` as a prefilter to find seeds, then computes a bounded
approximation of HMMER3 Forward/Backward around them.

## Branches

Small fixes go straight on `dev`. A feature or a big fix gets its own branch,
cut from `dev` and merged back into `dev` when it is done. A local merge is
fine for that.

`main` carries releases and nothing else. A release is a pull request from
`dev` to `main`. Never open that pull request, never merge it, and never push
to `main`. Releasing is not your job.

## Releases

Both crates follow [semantic versioning](https://semver.org/spec/v2.0.0.html).
While a crate is below 1.0, a breaking change takes the minor number.

The crates are versioned separately. Bump a crate only when it changed. When
libnail is bumped, change the `version` on the libnail dependency in
`nail/Cargo.toml` to match.

CI cuts the tags. Do not tag by hand. On a push to `main`,
`.github/workflows/release.yml` compares each crate's version with the previous
`main` commit, creates an annotated `nail-vX.Y.Z` or `libnail-vX.Y.Z` tag for
each crate whose version changed, and builds release binaries.

## Formatting

`rustfmt` with default settings is the format. Run `cargo fmt` before
committing.

## Changelog

Each crate has its own `CHANGELOG.md`, following
[Keep a Changelog](https://keepachangelog.com/en/1.0.0/). When a change alters
what a user of a crate sees, add an entry under `[Unreleased]` in that crate's
changelog, in the same commit. libnail's public API goes in
`libnail/CHANGELOG.md`. The CLI, its output, and nail's own modules go in
`nail/CHANGELOG.md`. A change that touches both crates gets an entry in each.

Entries match the ones already there:

- One bullet per public item, under `### Added`, `### Changed`, `### Removed`
  or `### Fixed`, in that order, with only the headings that have entries.
- Lowercase, no trailing period, opening with the past-tense verb of its
  section: `added`, `removed`, `fixed`. Under Changed, `renamed`, `moved`,
  `split`, `refactored`, or a plain statement such as
  "struct `Profile` now derives `PartialEq`".
- Name the kind of item, then its identifier in backticks:
  "added struct `AmbiguityMap`",
  "added field `use_accession` to struct `AlignConfig`",
  "added CLI params `--use-accession`, `--tbl-format`",
  "added mod `io::database`". Group members of one type with braces:
  "added methods `MmseqsDbPaths::{destroy(), check()}`".
- A refactor with several parts is one bullet with four-space-indented
  sub-bullets, one per part.
- In the commit that bumps the version, turn `[Unreleased]` into
  `## [X.Y.Z] - YYYY-M-D`, month and day unpadded, and leave a fresh empty
  `## [Unreleased]` above it.
- The `<!-- ************* -->` comment is a hand-placed marker. Do not move or
  remove it.

## Tests

`cargo test` runs unit tests; there are no integration tests. Fixtures live in
`fixtures/` at the workspace root, and tests reach them through
`CARGO_MANIFEST_DIR`. No test needs `mmseqs`. Both Cargo profiles build at
opt-level 3, so the first test build is slow.

## Platforms

Linux, macOS and Windows on x86_64 and aarch64; CI builds all of them. Nothing
in the source is target-specific. At runtime nail needs an MMseqs2 binary,
found through `--mmseqs-path` (default `mmseqs`).
