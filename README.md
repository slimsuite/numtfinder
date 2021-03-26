# NUMTFinder: Nuclear mitochondrial fragment (NUMT) search tool

```
NUMTFinder v0.3.0
```

For a better rendering and navigation of this document, please download and open [`./docs/numtfinder.docs.html`](./docs/numtfinder.docs.html), or visit <https://slimsuite.github.io/numtfinder/>.
Documentation can also be generated by running NUMTFinder with the `dochtml=T` option. (R and pandoc must be installed - see below.)

## Introduction

NUMTFinder uses a mitochondrial genome to search against genome assembly and identify putative NUMTs. NUMT fragments
are then combined into NUMT blocks based on proximity.

The general NUMTFinder workflow is:

1. Generate a double-copy linearised mtDNA sequence from the circular genome.
2. Perform a BLAST+ blastn search of the double-mtDNA versus the genome assembly using GABLAM.
3. Optionally filter short NUMT hits based on length.
4. Optionally filter NUMT hits based on hit sequence name and/or high identity (e.g. identify/remove real mtDNA).
5. Collapse nearby fragments into NUMT blocks. By default, fragments can incorporate duplications and rearrangements,
including inversions. Setting `stranded=T` will restrict blocks to fragments on the same strand.
6. Map fragments back on to the mtDNA genome and output a coverage plot.

Plans for future releases include:
* incorporation of additional search methods (LAST or kmers)
* assembly masking options
* options to restrict NUMT blocks to fully collinear hits.
* automated running of Diploidocus long-read regcheck on fragments and blocks

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

To generate mtDNA coverage plots, R will need to be installed. To generate documentation with `dochtml`, R will
need to be installed and a pandoc environment variable must be set, e.g.

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
exclude=LIST    : Exclude listed sequence names from search [mtDNA sequence name]
mtmaxcov=PERC   : Maximum percentage coverage of mtDNA (at mtmaxid identity) to allow [99]
mtmaxid=PERC    : Maximum percentage identity of mtDNA hits > mtmaxcov coverage to allow [99]
mtmaxexclude=T/F: Whether add sequences breaching mtmax filters to the exclude=LIST exclusion list [True]
keepblast=T/F   : Whether to keep the blast results files rather than delete them [True]
forks=INT       : Use multiple threads for the NUMT search [0]
### ~ NUMTFinder block options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
fragmerge=X     : Max Length of gaps between fragmented local hits to merge [8000]
stranded=T/F    : Whether to only merge fragments on the same strand [False]
### ~ NUMTFinder output options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
localgff=T/F    : Whether to output GFF format files of the NUMT hits against the genome [True]
localsam=T/F    : Whether to output SAM format files of the NUMT hits against the genome [True]
depthplot=T/F   : Whether to output mtDNA depth plots of sequence coverage (requires R) [True]
fasdir=PATH     : Directory in which to save fasta files [numtfasta/]
fragfas=T/F     : Whether to output NUMT fragment to fasta file [True]
fragrevcomp=T/F : Whether to reverse-complement DNA fragments that are on reverse strand to query [True]
blockfas=T/F    : Whether to generate a combined fasta file of NUMT block regions (positive strand) [True]
### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
```

---

# NUMTFinder Workflow

## NUMTFinder setup

NUMTFinder will start by checking and summarising the input sequences. NUMTFinder requires a fasta file containing
the genome to search for NUMTs (`seqin=FILE`) and the mitochondrial genome to search (`mtdna=FILE`), both in fasta
format. By default, NUMTFinder expects the mtDNA to be a single full-length circular sequence with no overhangs.
If this is not the case, switch `circle=F` and NUMTFinder will search all sequences in the file as simple linear
sequences. If either file is missing or fails to load a sequence, NUMTFinder will exit.

By default, any genome sequence matching the mtDNA sequence name(s) will be excluded from the search. This can be
over-ridden by setting the `exclude=LIST` option. Note that these sequences are still included in the BLAST
search and GABLAM summary files (including fragment SAM and GFF output), but will be excluded from the NUMT
fragment and block output, including the mtDNA coverage plot.

## mtDNA circularisation

If `circle=T` then NUMTFinder will generate a double-length mtDNA sequence for the actual search. This is to stop
artificial fragmentation of NUMTs across the circularisation breakpoint. A new double-length sequence will be
output to `$BASEFILE.mtdna2X.fasta` and used as the query for the [GABLAM](http://rest.slimsuite.unsw.edu.au/gablam)
search. This sequence will have "2X" appended to its sequence name. If `force=F` and this file already exists,
it will not be recreated.

## GABLAM (BLAST+) search

NUMTFinder uses [GABLAM](http://rest.slimsuite.unsw.edu.au/gablam) for the actual NUMT search. GABLAM is a wrapper
for [BLAST+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download) and, by default,
a `blastn` search is performed with an evalue threshold of `1e-4`. This can be changed with `blaste=X`. Other
BLAST+ settings can be modified as per [GABLAM](http://rest.slimsuite.unsw.edu.au/gablam) and [rje_blast](http://rest.slimsuite.unsw.edu.au/blast)
documentation. By default, hits of any length are accepted. This can be made more stringent with `minfraglen=INT`
to set the minimum fragment lenght in basepairs.

**NOTE:** The NUMTFinder defaults are set to be reasonably relaxed. In particular, it is possible that repeat
sequences in the mtDNA might result in multiple, quite short, hits in the search genome. This should be apparent
in the coverage plot (below).

GABLAM outputs a number of files with the prefix `$BASEFILE.numtsearch.*`, including the BLAST results file itself
(`$BASEFILE.numtsearch.blast`), unless `keepblast=F`. Fasta files of the hit fragments will also be output into
`numtfasta/` (`fasdir=PATH`). Parts of the GABLAM output can be switched off with `localsam=F`, `localgff=F` and/or
`fragfas=F`.


## NUMT filtering

Main processing of NUMTs by NUMTFinder uses the `$BASEFILE.numtsearch.unique.tdt` output from GABLAM. This reduces
the local BLAST hits to unique coverage of the `seqin=FILE` genome, by starting with the hits containing the
biggest number of identical matching bases and then progressively trimming overlapping hits until there are none.
This should reduce the NUMT fragments to the largest contiguous hits in the genome, although it is possible that
an unusual (and complex) mtDNA structure could result in expected behaviour and fragmentation of hits.

Many genome assemblies contain the mtDNA as well as the nuclear sequences. By default, NUMTFinder will look and
screen these out using quite stringent filters of >99% mtDNA coverage at >99% identity. Any sequences with 1+
hits at this stringency will be added to the `exclude=LIST` exclusion list, and reported. This can be relaxed to
only filter the hits themselves with `mtmaxexclude=F`. If a suspected mtDNA sequence contains more NUMT fragments,
a warning message will be added to the log.

Finally, any hits to sequence excluded by `exclude=LIST` will be removed from the results.

## NUMT fragments

Main NUMTFinder output is a table of predicted NUMT fragments, generated from the filtered GABLAM results. The
strand (`+`/`-`) is added and start/end positions modified to always be relative to the +ve strand. Unless
`circle=F`, any hits over-shooting the original circular mtDNA will have their positions modified to match the
original mtDNA sequence. For example, a match to region 16001-18000 of a double-length 17kb mtDNA will become a
match to region 16001-1000. Any matches where `mtEnd` is smaller than `mtStart` indicates a hit of this nature
that spans the circularisation point.

Fragments are then saved to `$BASEFILE.numtfrag.tdt`:

```
SeqName,Start,End,Strand,BitScore,Expect,Length,Identity,mtStart,mtEnd
```


## NUMT Blocks

Next, NUMTFinder merges nearby fragments into longer NUMT "blocks". This is done simply on proximity in the
assembly sequence. By default, NUMT fragments within 8kb of each other will be merged into a common block.
This is controlled by `fragmerge=INT` (bp). Where merged NUMT fragments are on both strands, `Strand` will be
set to `+/-`. The list of `FragNum` fragment positions will be output in the `mtFrag` field, with their combined
length in `FragLen`. Note that this is their length in the assembly, *not* the BLAST alignment `Length`, which is
summed in the `Length` field. The combined non-NUMT regions between merged fragments is given in `FragGaps`.

Blocks are then saved to `$BASEFILE.numtblock.tdt`:

```
SeqName,Start,End,Strand,BitScore,Expect,Length,Identity,mtFrag,FragNum,FragLen,FragGaps
```

**NOTE:** NUMT block fasta output is not yet implemented. Please contact the author if this is desired.

## Coverage plot

The final part of the NUMTFinder pipeline is to map all the (filtered) NUMT fragments back on to the original
mtDNA and generate plot of the depth of coverage by NUMT fragments across the mtDNA genome. Note that for this
output, any fragments spanning the circularisation point *will* be divided into two fragments. This will not
affect the depth plot, but will alter the accompanying "readlen" plot.


---

# NUMTFinder Outputs

NUMTFinder output will be named using the `basefile=X` prefix (hereon `$BASEFILE`), which is `numtfinder` by
default. The default NUMTFinder outputs are:

```
|-- numtfasta/
|   +-- $MTACC2X.fas
|-- $BASEFILE.coverage.tdt
|-- $BASEFILE.depthplot.tdt
|-- $BASEFILE.dirnlenplot.tdt
|-- $BASEFILE.log
|-- $BASEFILE.mtdna2X.acc.fas
|-- $BASEFILE.mtdna2X.fasta
|-- $BASEFILE.mtdna2X.fasta.index
|-- $BASEFILE.numtblock.tdt
|-- $BASEFILE.numtfrag.tdt
|-- $BASEFILE.numtsearch.blast
|-- $BASEFILE.numtsearch.gablam.tdt
|-- $BASEFILE.numtsearch.hitsum.tdt
|-- $BASEFILE.numtsearch.local.gff
|-- $BASEFILE.numtsearch.local.sam
|-- $BASEFILE.numtsearch.local.tdt
|-- $BASEFILE.numtsearch.unique.gff
|-- $BASEFILE.numtsearch.unique.sam
|-- $BASEFILE.numtsearch.unique.tdt
|-- $BASEFILE.readlenplot.tdt
|-- $BASEFILE.readlen.tdt
|-- $BASEFILE.rid.tdt
+-- $BASEFILE.SAMPlots/
    |-- $BASEFILE.depth.$MTACC.png
    +-- $BASEFILE.readlen.$MTACC.png
```

The main NUMTFinder outputs are:

* $BASEFILE.log
* $BASEFILE.numtsearch.unique.gff
* $BASEFILE.numtsearch.unique.sam
* $BASEFILE.numtblock.tdt
* $BASEFILE.numtfrag.tdt
* $BASEFILE.SAMPlots/$BASEFILE.depth.$MTACC.png
* numtfasta/$MTACC2X.fas

More details will be added in future releases.



<br>
<small>&copy; 2021 Richard Edwards | richard.edwards@unsw.edu.au</small>
