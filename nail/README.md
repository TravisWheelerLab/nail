# nail

## About

nail is a biological sequence alignment tool. Currently, only protein sequence alignment is supported.

nail uses [MMseqs2](https://github.com/soedinglab/MMseqs2) to find rough alignment seeds to use as a starting point
for a highly sensitive, bounded sequence alignment algorithm.

## Installation

To build nail from source, you'll first need to install Rust and Cargo.
The easiest way to do that is to use [rustup](https://rustup.rs/).

Once that's done, you can then build nail:

    git clone https://github.com/TravisWheelerLab/nail
    cd nail/
    cargo build --release

You'll then find the compiled binary at: `target/release/nail`

For example, try running:

    target/release/nail -h

## Usage

The input to nail is a query multiple sequence alignment (stockholm) file and a target sequence (fasta) file.
To run the nail pipeline, use the `nail search` command:

For example:

    $ nail search query.sto target.fa

## License

nail is licensed under the BSD-3-Clause license.

See `LICENSE` for details.

## Authors

Jack Roddy - jroddy@arizona.edu
