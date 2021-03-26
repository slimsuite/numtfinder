# NUMTFinder: Nuclear mitochondrial fragment (NUMT) search tool

```
NUMTFinder v0.2.0
```

For a better rendering and navigation of this document, please download and open [`./docs/numtfinder.docs.html`](./docs/numtfinder.docs.html), or visit <https://slimsuite.github.io/numtfinder/>.
Documentation can also be generated by running NUMTFinder with the `dochtml=T` option. (R and pandoc must be installed - see below.)

## Introduction

NUMTFinder uses a mitochondrial genome to search against genome assembly and identify putative NUMTs. NUMT fragments
are then combined into NUMT blocks based on proximity.

The general NUMTFinder workflow is:

1. Generate a double-copy linearised mtDNA sequence from the circular genome.
2. Perform a BLAST+ blastn search of the double-mtDNA versus the genome assembly using GABLAM.
3. Optionally filter NUMT hits based on length.
4. Collapse nearby fragments into NUMT blocks. By default, fragments can incorporate duplications and rearrangements,
including inversions. Setting `stranded=T` will restrict blocks to fragments on the same strand.

Plans for future releases include:
* incorporation of additional search methods (LAST or kmers)
* assembly masking options
* options to restrict NUMT blocks to fully collinear hits.
* automated running of Diploidocus regcheck on fragments and blocks
* depth profile of coverage across mitochondrion

---

# Running NUMTFinder

NUMTFinder is written in Python 2.x and can be run directly from the commandline:

    python $CODEPATH/numtfinder.py [OPTIONS]

If running as part of [SLiMSuite](http://slimsuite.blogspot.com/), `$CODEPATH` will be the SLiMSuite `tools/`
directory. If running from the standalone [NUMTFinder git repo](https://github.com/slimsuite/numtfinder), `$CODEPATH`
will be the path the to `code/` directory. Please see details in the [NUMTFinder git repo](https://github.com/slimsuite/numtfinder)
for running on example data.

## Dependencies

BLAST+ must be installed and either added to the environment `$PATH` or given to NUMTFinder with the `blast+path` setting.

To generate documentation with `dochtml`, R will need to be installed and a pandoc environment variable must be set, e.g.

    export RSTUDIO_PANDOC=/Applications/RStudio.app/Contents/MacOS/pandoc

For NUMTFinder documentation, run with `dochtml=T` and read the `*.docs.html` file generated.

## Commandline options

A list of commandline options can be generated at run-time using the `-h` or `help` flags. Please see the general
[SLiMSuite documentation](http://slimsuite.blogspot.com/2013/08/command-line-options.html) for details of how to
use commandline options, including setting default values with **INI files**.

```
### ~ Main NUMTFinder run options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
seqin=FILE      : Genome assembly in which to search for NUMTs []
mtdna=FILE      : mtDNA reference genome to use for search []
basefile=X      : Prefix for output files [numtfinder]
summarise=T/F   : Whether to summarise input sequence files upon loading [True]
dochtml=T/F     : Generate HTML NUMTFinder documentation (*.docs.html) instead of main run [False]
### ~ NUMTFinder search options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
circle=T/F      : Whether the mtDNA is circular [True]
blaste=X        : BLAST+ blastn evalue cutoff for NUMT search [1e-4]
minfraglen=INT  : Minimum local (NUMT fragment) alignment length (sets GABLAM localmin=X) [0]
keepblast=T/F   : Whether to keep the blast results files rather than delete them [True]
forks=INT       : Use multiple threads for the NUMT search [0]
### ~ NUMTFinder block options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
fragmerge=X     : Max Length of gaps between fragmented local hits to merge [8000]
stranded=T/F    : Whether to only merge fragments on the same strand [False]
### ~ Sequence output options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
fasdir=PATH     : Directory in which to save fasta files [numtfasta/]
fragfas=T/F     : Whether to output NUMT fragment to fasta file [False]
fragrevcomp=T/F : Whether to reverse-complement DNA fragments that are on reverse strand to query [True]
blockfas=T/F    : Whether to generate a combined fasta file of NUMT block regions (positive strand) [True]
### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
```


<br>
<small>&copy; 2021 Richard Edwards | richard.edwards@unsw.edu.au</small>
