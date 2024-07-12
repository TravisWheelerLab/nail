# nail workspace

## About

This is a cargo workspace for nail, which is a profile Hidden Markov Model (pHMM) biological sequence alignment tool.
Using the fast [MMseqs2](https://github.com/soedinglab/MMseqs2) search pipeline to produce candidate alignment seeds, nail computes a fast approximation of the [HMMER3](http://hmmer.org/) Forward/Backward (F/B) sequence alignment algorithm.
Currently, nail only supports amino acid search, with nucleotide search coming in a future update. 

### What's here
There are two sub-projects in the nail workspace:
1. `nail`: this is the command line tool
2. `libnail`: this is a Rust library that contains the implementation of nail's sparse alignment algorithms

## Dependencies

The `nail search` pipeline uses the `mmseqs search` tool as an alignment prefilter.
In the future, we will replace this with our own prefiltering strategies.

nail has been tested with [MMseqs2 Release 15-6f452](https://github.com/soedinglab/MMseqs2/releases/tag/15-6f452).
We have not tested nail against other versions of MMseqs2, but they may work.

To run the `nail search` pipeline, `mmseqs search` must be available in your system path.

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

The nail command line interface has two subcommands: `search` and `seed`.

### nail search

The `nail search` command runs the entire nail pipeline, including running MMseqs2 to find alignment seeds.

The input to `nail search` is a query file (.hmm|.fasta) and a target sequence database file (.fasta).

By default, the search results will be written to `./results.tbl` in a tabular format, and alignment output is written to stdout.
In addition, a collection of temporary files required to run `mmseqs search`, will be written to the `./prep/` directory.

For example, when running nail search, you'll see something like:

    $ nail search query.hmm target.fa

    ==  score: 257.1 bits;  E-value: 9.7e-79
                                  7tm_1     2 NllVilvilrnkklrtptnifllnLavaDllvlllvlpfslvyallegdwvfgevlCklvtaldvvnltasillltaisi 81   
                                              N +Vi++++r++kl+tp+n+++ +La++Dllv++lv+p+s++y++ +g+w++g++lC+++ + d++++tasi++l++i++
    O08892|reviewed|5-hydroxytryptamine    66 NAFVIATVYRTRKLHTPANYLIASLAFTDLLVSILVMPISTMYTV-TGRWTLGQALCDFWLSSDITCCTASIMHLCVIAL 145  
                                              9********************************************06*********************************
    
                                  7tm_1    82 DRYlaIvkplkykrirtkrralvlilvvWvlalllslppllfsgtktesaekeetvClidfpeeestwevsytlllsvlg 161  
                                              DRY+aI+ ++ y+++rt+rra+ +i++vWv+++++slpp+++++ k e   +e+  Cl+++++      v yt++++ ++
    O08892|reviewed|5-hydroxytryptamine   146 DRYWAITDAVGYSAKRTPRRAAGMIALVWVFSICISLPPFFWRQAKAE---EEVLDCLVNTDH------VLYTVYSTGGA 225  
                                              *****************************************777666500099******99990000009**********
    
                                  7tm_1   162 fllpllvilvcyvrilrtlrksakkeks.................................................... 241  
                                              f+lp+l+++ +y ri+ ++r++  k++                                                     
    O08892|reviewed|5-hydroxytryptamine   226 FYLPTLLLIALYGRIYVEARSRILKQTPNKTGKRLTRAQLITDSPGSTSSVTSINSRAPEVPCDSGSPVYVNQVKVRVSD 305  
                                              ***********************9999899999999999999999999********************************
    
                                  7tm_1   242 ....rkkksarkerkalktllvvvvvfvlcwlPyfilllldsllkeceseklve.tallitlllayvnsclNPiiY 317  
                                                  +kk +a++erka+ktl+v++++f++cwlP+fi++l++ +   c++ +  + ++++++++l+y+ns++NPiiY
    O08892|reviewed|5-hydroxytryptamine   306 ALLEKKKLMAARERKATKTLGVILGAFIVCWLPFFIISLVMPI---CKDACWFHMAIFDFFTWLGYLNSLINPIIY 381  
                                          ***9999************************************000777666551666******************
    ...
    ...
    ...

Checking the `results.tbl` might look something like:

    $ head -n 15 results.bl

    #                                                                         target target query query       comp         cell  
    # target                                                  query           start  end    start end   score bias evalue  frac  
    # ------------------------------------------------------- --------------- ------ ------ ----- ----- ----- ---- ------- ----- 
    F1MV99|reviewed|Somatostatin                              7tm_1-consensus 58     306    1     260   175.8 7.5  1.0e-53 0.066
    O08858|reviewed|Somatostatin                              7tm_1-consensus 54     303    1     260   173.0 7.8  7.7e-53 0.069
    C3ZQF9|reviewed|QRFP-like                                 7tm_1-consensus 64     326    1     260   149.3 8.6  1.3e-45 0.070
    O02813|reviewed|Neuropeptide                              7tm_1-consensus 56     319    1     260   147.9 5.4  3.5e-45 0.070
    O08565|reviewed|C-X-C                                     7tm_1-consensus 52     299    1     260   145.8 5.3  1.5e-44 0.076
    O02835|reviewed|Neuropeptide                              7tm_1-consensus 57     320    1     260   145.7 6.8  1.6e-44 0.071
    A0T2N3|reviewed|Apelin                                    7tm_1-consensus 51     316    1     260   145.0 4.3  2.6e-44 0.077
    O08726|reviewed|Galanin                                   7tm_1-consensus 42     292    1     260   144.9 4.2  2.8e-44 0.073
    O08556|reviewed|C-C                                       7tm_1-consensus 49     299    1     260   144.3 5.3  4.4e-44 0.073
    F1R332|reviewed|Galanin                                   7tm_1-consensus 35     286    1     260   143.3 3.8  8.9e-44 0.077
    D4A7K7|reviewed|G-protein                                 7tm_1-consensus 44     304    1     260   142.7 4.3  1.3e-43 0.073
    E7F7V7|reviewed|Galanin                                   7tm_1-consensus 43     294    1     260   142.6 1.7  1.4e-43 0.073


### nail seed

The `nail seed` command runs MMseqs2 and produces a `seeds.json` file.

For example:

    $ nail seed query.hmm target.fa

Seeds can be provided to `nail search` using the `--seeds <seeds.json>` flag, which will skip the seed step in the search pipeline.

    $ nail search --seeds seeds.json query.hmm target.fa

In practice, these seeds may be produced from any source as long as they are formatted in the following way:

```
{
  "query1": {
    "target1": {
      "target_start": 48,
      "target_end": 287,
      "profile_start": 1,
      "profile_end": 259,
      "score": 168.0
    },
    "target2": {
      "target_start": 72,
      "target_end": 343,
      "profile_start": 23,
      "profile_end": 259,
      "score": 106.0
    },
  "query2": {
    "target3": {
      "target_start": 56,
      "target_end": 303,
      "profile_start": 1,
      "profile_end": 259,
      "score": 125.0
    },
  }
  ...
}
```

## License

nail is licensed under the BSD-3-Clause license.

See `LICENSE` for details.

## Authors

Jack Roddy - jroddy@arizona.edu
