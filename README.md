# nail workspace

## About

This is a cargo workspace for nail, which is a profile Hidden Markov Model (pHMM) biological sequence alignment tool.
Using the fast [MMseqs2](https://github.com/soedinglab/MMseqs2) search pipeline to produce candidate alignment seeds, nail computes a fast approximation of the [HMMER3](http://hmmer.org/) Forward/Backward (F/B) sequence alignment algorithm.
Currently, nail only supports amino acid search, with nucleotide search coming in a future update.

The nail preprint paper can be found on bioRxiv (doi: https://doi.org/10.1101/2024.01.27.577580).

### What's here
There are two sub-projects in the nail workspace:
1. `nail`: this is the command line tool
2. `libnail`: this is a Rust library that contains the implementation of nail's sparse alignment algorithms

## Dependencies

The `nail search` pipeline uses `mmseqs search` from the MMseqs2 suite and the `hmmbuild` tool from the HMMER3 suite.
In the future, we plan to replace both dependencies with our own alternatives that will be internal to nail.

nail has been tested with [MMseqs2 Release 15-6f452](https://github.com/soedinglab/MMseqs2/releases/tag/15-6f452) and [HMMER v3.4](http://hmmer.org/download.html).
We have not tested nail against versions of these tools, but they may work.

To run the `nail search` pipeline, `mmseqs search` and `hmmbuild` must be available in your system path.

## Installation

To install nail, you'll need the Rust development tooling, which means you'll need to install the Rust compiler and Cargo.
The easiest way to do that is to use [rustup](https://rustup.rs/).

Once Cargo is installed, you can install nail with:

    cargo install nail

### Building from source

If you'd like to build nail from source, you can clone this repository and build the project with Cargo:

    git clone https://github.com/TravisWheelerLab/nail
    cd nail/
    cargo build --release

You'll then find the compiled binary at: `target/release/nail`

For example, try running:

    target/release/nail -h

## Usage

The nail command line interface has four subcommands: `search`, `prep`, `seed`, and `align`.

### Running the nail pipeline with nail search

`nail search` runs each step of the nail pipeline in succession: `prep`, `seed`, and `align`.

The input to `nail search` is a query multiple sequence alignment (MSA) file in the Stockholm format and a target sequence database file in the FASTA format.

For example, running

    $ nail search query.sto target.fa

By default, the search results will be written to `./results.tsv` in a tabular format, and alignment output is written to stdout.
In addition, a collection of temporary files required to run `mmseqs search`, along with a pHMM built from `query.sto`, will be written to the `./prep/` directory.

### Running individual pipeline steps

nail also supports running each pipeline step in isolation.
This may be particularly useful if:
- you'd like to cache the results of either `prep` or `seed` and experiment with parameters on a successive step, or
- you'd like to experiment with alignment seeds produced by a source other than MMseqs2

#### nail prep

Running

    $ nail prep query.sto target.fa

will create a `./prep` directory including the files that are used in `nail seed`.

The prep directory is full of the *many* custom format input files used by the `MMseqs2 search` tool and the query pHMM created with HMMER3's `hmmbuild` tool.

For example:

```
$ tree prep/
prep/
├── msaDB
├── msaDB.dbtype
├── msaDB.index
├── query.hmm
├── queryDB
├── queryDB.dbtype
├── queryDB.index
├── queryDB_h
├── queryDB_h.dbtype
├── queryDB_h.index
├── targetDB
├── targetDB.dbtype
├── targetDB.index
├── targetDB.lookup
├── targetDB.source
├── targetDB_h
├── targetDB_h.dbtype
└── targetDB_h.index
```
#### nail seed

Running 

    $ nail seed prep/

where `prep/` is a directory created as result of running `nail prep`, will produce the file `./seeds.json`, which has the following structure:

```
{
  "7tm_1": [
    {
      "target_name": "P07700|reviewed|Beta-1",
      "target_start": 58,
      "target_end": 343,
      "profile_start": 1,
      "profile_end": 260
    },
    {
      "target_name": "P34971|reviewed|Beta-1",
      "target_start": 75,
      "target_end": 366,
      "profile_start": 1,
      "profile_end": 260
    },
    ...
}
```

#### nail align

Running 

    $ nail align query.hmm target.fa seeds.json

will run nail's sparse Forward/Backward alignment algorithm, producing optimal-scoring alignments.
By default, alignment output will be written to stdout, e.g:

```
==  score: 255.1 bits;  E-value: 7.5e-77
                 7tm_1     1 gNllVilvilrnkklrtptnifllnLavaDllvlllvlpfslvyallegdwvfgevlCklvtaldvvnltasillltais 80   
                             gN+lVi++i r+++l+t tn+f+++La+aDl+++llv+pf ++  + +g+w +g++lC+++t+ldv+++tasi +l++i+
P07700|reviewed|Beta-1    58 GNVLVIAAIGRTQRLQTLTNLFITSLACADLVMGLLVVPFGATLVV-RGTWLWGSFLCECWTSLDVLCVTASIETLCVIA 137  
                             0000000000000000000000000000000000*************************08*******************

                 7tm_1    81 iDRYlaIvkplkykrirtkrralvlilvvWvlalllslppllfsgtktesae.....keetvClidfpeeestwevsytl 160  
                             iDRYlaI+ p++y++++t+ ra+v+i+ vW++++l+s++p+++ ++++e+++     ++   C++ ++       + y++
P07700|reviewed|Beta-1   138 IDRYLAITSPFRYQSLMTRARAKVIICTVWAISALVSFLPIMMHWWRDEDPQALKCYQDPGCCDFVTN-------RAYAI 217  
                             **********************776608****************************************************

                 7tm_1   161 llsvlgfllpllvilvcyvrilrtlrksakkeks.................................rkkksarkerkal 240  
                              +s+++f++pll+++++y r++r+++++ +k ++                                 +++ +a +e+kal
P07700|reviewed|Beta-1   218 ASSIISFYIPLLIMIFVYLRVYREAKEQIRKIDRCEGRFYGSQEQPQPPPLPQHQPILGNGRASKRKTSRVMAMREHKAL 297  
                             ***************************99866666655555555555500000005************************

                 7tm_1   241 ktllvvvvvfvlcwlPyfilllldsllkeceseklvetallitlllayvnsclNPii 297  
                             ktl+++++vf+lcwlP+f++++++++    +++ ++++++ ++++l+y+ns++NPii
P07700|reviewed|Beta-1   298 KTLGIIMGVFTLCWLPFFLVNIVNVF----NRDLVPDWLFVFFNWLGYANSAFNPII 354  
                             *************9*************************999998777777788888
```

By default, tabular results will be placed in `results.tsv`.
The default column order is:

1. target name
2. query name
3. target start
4. target end
5. query start
6. query end
7. bit score
8. composition bias score
9. E-value
10. cell fraction (the fraction of cells computed by sparse F/B)

For example:
```
$ head -n 5 results.tsv

P07700|reviewed|Beta-1 7tm_1 58 342 1 259 255.08 17.60 7.5e-77 6.4e-2
P34971|reviewed|Beta-1 7tm_1 75 365 1 259 254.19 12.99 1.4e-76 7.1e-2
P18090|reviewed|Beta-1 7tm_1 75 365 1 259 251.74 14.10 7.8e-76 7.1e-2
Q28927|reviewed|Beta-1 7tm_1 75 363 1 259 250.06 16.87 2.5e-75 7.0e-2
```

In principle, the seeds used by `nail align` could come from any source.
If you'd like to test out your own prefilter & seeding methods using nail's sparse F/B alignment stage, you simply need to produce a `seeds.json` file that conforms the structure described above.

## License

nail is licensed under the BSD-3-Clause license.

See `LICENSE` for details.

## Authors

Jack Roddy - jroddy@arizona.edu
