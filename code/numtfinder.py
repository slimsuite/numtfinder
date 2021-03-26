#!/usr/bin/python

# See below for name and description
# Copyright (C) 2016 Richard J. Edwards <redwards@cabbagesofdoom.co.uk>
#  
# This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
# warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with this program; if not, write to 
# the Free Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
#
# Author contact: <seqsuite@gmail.com> / School of Biotechnology and Biomolecular Sciences, UNSW, Sydney, Australia.
#
# To incorporate this module into your own programs, please see GNU Lesser General Public License disclaimer in rje.py

"""
Module:       NUMTFinder
Description:  Nuclear mitochondrial fragment (NUMT) search tool
Version:      0.3.0
Last Edit:    26/03/21
Citation:     Edwards RJ et al. (2021), BMC Genomics [PMID: 33726677]
GitHub:       https://github.com/slimsuite/numtfinder
Copyright (C) 2021  Richard J. Edwards - See source code for GNU License Notice

Function:
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

Commandline:
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
"""
#########################################################################################################################
### SECTION I: GENERAL SETUP & PROGRAM DETAILS                                                                          #
#########################################################################################################################
import os, string, sys, time
slimsuitepath = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../')) + os.path.sep
sys.path.append(os.path.join(slimsuitepath,'libraries/'))
sys.path.append(os.path.join(slimsuitepath,'tools/'))
### User modules - remember to add *.__doc__ to cmdHelp() below ###
import rje, rje_db, rje_obj, rje_rmd, rje_seqlist, rje_samtools
import gablam
#########################################################################################################################
def history():  ### Program History - only a method for PythonWin collapsing! ###
    '''
    # 0.0.0 - Initial Compilation.
    # 0.1.0 - Added dochtml=T and modified docstring for standalone git repo.
    # 0.1.1 - Fixed bug with default fragmerge=INT. Now set to 8kb.
    # 0.2.0 - Added SAM output and depth profile of coverage across mitochondrion.
    # 0.3.0 - Added additional exclusion, flagging and filtering of possible mtDNA.
    '''
#########################################################################################################################
def todo():     ### Major Functionality to Add - only a method for PythonWin collapsing! ###
    '''
    # [Y] : Populate Module Docstring with basic info.
    # [Y] : Populate makeInfo() method with basic info.
    # [Y] : Add full description of program to module docstring.
    # [?] : Create initial working version of program.
    # [ ] : Add REST outputs to restSetup() and restOutputOrder()
    # [ ] : Add to SLiMSuite or SeqSuite.
    # [ ] : Add function to give other genomes to tile and then analyse for coverage with Diploidocus.
    # [ ] : Check and replace forks=INT with threads=INT.
    # [ ] : Add long read coverage analysis.
    # [Y] : Depth profile of coverage across mitochondrion
    # [ ] : Add option to exclude certain sequences (e.g. the mitochondrion!)
    # [ ] : Add option to identify and exclude hits above a certain length and identity (e.g. the mitochondrion!)
    '''
#########################################################################################################################
def makeInfo(): ### Makes Info object which stores program details, mainly for initial print to screen.
    '''Makes Info object which stores program details, mainly for initial print to screen.'''
    (program, version, last_edit, copy_right) = ('NUMTFinder', '0.3.0', 'March 2021', '2021')
    description = 'Nuclear mitochondrial fragment (NUMT) search tool'
    author = 'Dr Richard J. Edwards.'
    comments = ['This program is still in development and has not been published.',rje_obj.zen()]
    return rje.Info(program,version,last_edit,description,author,time.time(),copy_right,comments)
#########################################################################################################################
def cmdHelp(info=None,out=None,cmd_list=[]):   ### Prints *.__doc__ and asks for more sys.argv commands
    '''Prints *.__doc__ and asks for more sys.argv commands.'''
    try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        if not info: info = makeInfo()
        if not out: out = rje.Out()
        ### ~ [2] ~ Look for help commands and print options if found ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        cmd_help = cmd_list.count('help') + cmd_list.count('-help') + cmd_list.count('-h')
        if cmd_help > 0:
            rje.printf('\n\nHelp for {0} {1}: {2}\n'.format(info.program, info.version, time.asctime(time.localtime(info.start_time))))
            out.verbose(-1,4,text=__doc__)
            if rje.yesNo('Show general commandline options?',default='N'): out.verbose(-1,4,text=rje.__doc__)
            if rje.yesNo('Quit?'): sys.exit()           # Option to quit after help
            cmd_list += rje.inputCmds(out,cmd_list)     # Add extra commands interactively.
        elif out.stat['Interactive'] > 1: cmd_list += rje.inputCmds(out,cmd_list)    # Ask for more commands
        ### ~ [3] ~ Return commands ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        return cmd_list
    except SystemExit: sys.exit()
    except KeyboardInterrupt: sys.exit()
    except: rje.printf('Major Problem with cmdHelp()')
#########################################################################################################################
def setupProgram(): ### Basic Setup of Program when called from commandline.
    '''
    Basic Setup of Program when called from commandline:
    - Reads sys.argv and augments if appropriate
    - Makes Info, Out and Log objects
    - Returns [info,out,log,cmd_list]
    '''
    try:### ~ [1] ~ Initial Command Setup & Info ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        info = makeInfo()                                   # Sets up Info object with program details
        if len(sys.argv) == 2 and sys.argv[1] in ['version','-version','--version']: rje.printf(info.version); sys.exit(0)
        if len(sys.argv) == 2 and sys.argv[1] in ['details','-details','--details']: rje.printf('%s v%s' % (info.program,info.version)); sys.exit(0)
        if len(sys.argv) == 2 and sys.argv[1] in ['description','-description','--description']: rje.printf('%s: %s' % (info.program,info.description)); sys.exit(0)
        cmd_list = rje.getCmdList(sys.argv[1:],info=info)   # Reads arguments and load defaults from program.ini
        out = rje.Out(cmd_list=cmd_list)                    # Sets up Out object for controlling output to screen
        out.verbose(2,2,cmd_list,1)                         # Prints full commandlist if verbosity >= 2 
        out.printIntro(info)                                # Prints intro text using details from Info object
        cmd_list = cmdHelp(info,out,cmd_list)               # Shows commands (help) and/or adds commands from user
        log = rje.setLog(info,out,cmd_list)                 # Sets up Log object for controlling log file output
        return (info,out,log,cmd_list)                      # Returns objects for use in program
    except SystemExit: sys.exit()
    except KeyboardInterrupt: sys.exit()
    except: rje.printf('Problem during initial setup.'); raise
#########################################################################################################################
### END OF SECTION I                                                                                                    #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION II: New Class                                                                                               #
#########################################################################################################################
class NUMTFinder(rje_obj.RJE_Object):
    '''
    NUMTFinder Class. Author: Rich Edwards (2021).

    Str:str
    - mtDNA=FILE      : mtDNA reference genome to use for search []
    - mtQuery         : Actual mtDNA fasta file to use for search (double length if circle=T)
    - FasDir=PATH     : Directory in which to save fasta files [numtfasta/]
    - SeqIn=FILE      : Genome assembly in which to search for NUMTs []

    Bool:boolean
    - BlockFas=T/F    : Whether to generate a combined fasta file of NUMT block regions (positive strand) [True]
    - Circle=T/F      : Whether the mtDNA is circular [True]
    - DocHTML=T/F     : Generate HTML Snapper documentation (*.info.html) instead of main run [False]
    - FragFas=T/F     : Whether to output NUMT fragment to fasta file [False]
    - FragRevComp=T/F : Whether to reverse-complement DNA fragments that are on reverse strand to query [True]
    - MTMaxExclude=T/F: Whether add sequences breaching mtmax filters to the exclude=LIST exclusion list [True]
    - Stranded=T/F    : Whether to only merge fragments on the same strand [False]

    Int:integer
    - FragMerge=INT   : Max Length of gaps between fragmented local hits to merge [8000]
    - MinFragLen=INT  : Minimum local (NUMT fragment) alignment length (sets GABLAM localmin=X) [0]

    Num:float
    - MTMaxCov=PERC   : Maximum percentage coverage of mtDNA (at mtmaxid identity) to allow [99]
    - MTMaxID=PERC    : Maximum percentage identity of mtDNA hits > mtmaxcov coverage to allow [99]


    File:file handles with matching str filenames
    
    List:list
    - Exclude=LIST    : Exclude listed sequence names from search [mtDNA sequence name]

    Dict:dictionary    

    Obj:RJE_Objects
    - DB: Database object
    - GABLAM: GABLAM object
    - mtDNA: SeqList mtDNA reference genome to use for search []
    - SeqIn: SeqList genome assembly in which to search for NUMTs []
    '''
#########################################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes                                                                  #
#########################################################################################################################
    def _setAttributes(self):   ### Sets Attributes of Object
        '''Sets Attributes of Object.'''
        ### ~ Basics ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self.strlist = ['mtDNA','mtQuery','FasDir','SeqIn']
        self.boollist = ['BlockFas','Circle','DocHTML','FragFas','FragRevComp','MTMaxExclude','Stranded']
        self.intlist = ['FragMerge','MinFragLen']
        self.numlist = ['MTMaxCov','MTMaxID']
        self.filelist = []
        self.listlist = ['Exclude']
        self.dictlist = []
        self.objlist = ['DB','GABLAM','mtDNA','SeqIn']
        ### ~ Defaults ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self._setDefaults(str='None',bool=False,int=0,num=0.0,obj=None,setlist=True,setdict=True,setfile=True)
        self.setStr({})
        self.setBool({'BlockFas':True,'Circle':True,'DocHTML':False,'FragFas':False,'FragRevComp':True,'MTMaxExclude':True,'Stranded':False})
        self.setInt({'FragMerge':8000})
        self.setNum({'MTMaxCov':99.0,'MTMaxID':99.0})
        self.list['Exclude'] = 'mtDNA'
        ### ~ Other Attributes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self._setForkAttributes()   # Delete if no forking
        self.obj['DB'] = rje_db.Database(self.log,self.cmd_list+['tuplekeys=T'])
#########################################################################################################################
    def _cmdList(self):     ### Sets Attributes from commandline
        '''
        Sets attributes according to commandline parameters:
        - see .__doc__ or run with 'help' option
        '''
        for cmd in self.cmd_list:
            try:
                self._generalCmd(cmd)   ### General Options ### 
                self._forkCmd(cmd)  # Delete if no forking
                ### Class Options (No need for arg if arg = att.lower()) ### 
                #self._cmdRead(cmd,type='str',att='Att',arg='Cmd')  # No need for arg if arg = att.lower()
                #self._cmdReadList(cmd,'str',['Att'])   # Normal strings
                self._cmdReadList(cmd,'path',['FasDir'])  # String representing directory path
                self._cmdReadList(cmd,'file',['mtDNA','SeqIn'])  # String representing file path
                #self._cmdReadList(cmd,'date',['Att'])  # String representing date YYYY-MM-DD
                self._cmdReadList(cmd,'bool',['BlockFas','Circle','DocHTML','FragFas','FragRevComp','MTMaxExclude','Stranded'])  # True/False Booleans
                self._cmdReadList(cmd,'int',['FragMerge','MinFragLen'])   # Integers
                #self._cmdReadList(cmd,'float',['Att']) # Floats
                self._cmdReadList(cmd,'per',['MTMaxCov','MTMaxID']) # Percentages, stored 0-100 (<1 = x100)
                #self._cmdReadList(cmd,'min',['Att'])   # Integer value part of min,max command
                #self._cmdReadList(cmd,'max',['Att'])   # Integer value part of min,max command
                self._cmdReadList(cmd,'list',['Exclude'])  # List of strings (split on commas or file lines)
                #self._cmdReadList(cmd,'clist',['Att']) # Comma separated list as a *string* (self.str)
                #self._cmdReadList(cmd,'glist',['Att']) # List of files using wildcards and glob
                #self._cmdReadList(cmd,'cdict',['Att']) # Splits comma separated X:Y pairs into dictionary
                #self._cmdReadList(cmd,'cdictlist',['Att']) # As cdict but also enters keys into list
            except: self.errorLog('Problem with cmd:%s' % cmd)
#########################################################################################################################
    ### <2> ### Main Class Backbone                                                                                     #
#########################################################################################################################
    def run(self):  ### Main run method
        '''
        # NUMTFinder: Nuclear mitochondrial fragment (NUMT) search tool


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

        ```SeqName,Start,End,Strand,BitScore,Expect,Length,Identity,mtStart,mtEnd```


        ## NUMT Blocks

        Next, NUMTFinder merges nearby fragments into longer NUMT "blocks". This is done simply on proximity in the
        assembly sequence. By default, NUMT fragments within 8kb of each other will be merged into a common block.
        This is controlled by `fragmerge=INT` (bp). Where merged NUMT fragments are on both strands, `Strand` will be
        set to `+/-`. The list of `FragNum` fragment positions will be output in the `mtFrag` field, with their combined
        length in `FragLen`. Note that this is their length in the assembly, *not* the BLAST alignment `Length`, which is
        summed in the `Length` field. The combined non-NUMT regions between merged fragments is given in `FragGaps`.

        Blocks are then saved to `$BASEFILE.numtblock.tdt`:

        ```SeqName,Start,End,Strand,BitScore,Expect,Length,Identity,mtFrag,FragNum,FragLen,FragGaps```

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

        '''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.getBool('DocHTML'): return rje_rmd.docHTML(self)
            #i# Load the mtDNA and reference genome
            if not self.setup(): return False
            ### ~ [2] ~ Add main run code here ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [2a] ~ If circular, generate double-length mtDNA sequence ~~~~~~~~~~~~~~~~~~~~~~~~ ##
            self.mtQuery()
            ## ~ [2b] ~ Perform NUMT search (GABLAM by default) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            self.numtSearch()
            ## ~ [2c] ~ Load in NUMT results as numtfrag and merge to numtblock ~~~~~~~~~~~~~~~~~~~ ##
            self.numtProcess()
            ## ~ [2d] ~ Output NUMT Block sequences ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if self.getBool('BlockFas'):
                self.printLog('#SORRY','NUMT block fasta output not yet implemented!')
            ## ~ [2e] ~ Coverage output ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            self.coverageOutputs()
            return
        except:
            self.errorLog(self.zen())
            raise   # Delete this if method error not terrible
#########################################################################################################################
    def setup(self):    ### Main class setup method.
        '''Main class setup method.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.obj['mtDNA'] = rje_seqlist.SeqList(self.log,['summarise=T']+self.cmd_list+['seqin={0}'.format(self.getStr('mtDNA')),'dna'])
            if not self.obj['mtDNA'].seqNum(): raise IOError('Failed to load sequences from mtdna=FILE')
            elif self.obj['mtDNA'].seqNum() > 1 and self.getBool('Circle'):
                self.warnLog('{0} sequences loaded from mtdna=FILE and circle=T: will use first sequence only'.format(self.obj['mtDNA'].seqNum()))
            if self.list['Exclude'] == 'mtDNA':
                self.list['Exclude'] = self.obj['mtDNA'].names()
                if self.obj['mtDNA'].seqNum() > 1:
                    self.printLog('#MTDNA','Set mtDNA exclusion sequence names: {0}'.format('; '.join(self.list['Exclude'])))
                else:
                    self.printLog('#MTDNA','Set mtDNA exclusion sequence name: {0}'.format('; '.join(self.list['Exclude'])))
            self.obj['SeqIn'] = rje_seqlist.SeqList(self.log,['summarise=T']+self.cmd_list+['seqin={0}'.format(self.getStr('SeqIn')),'dna'])
            if not self.obj['SeqIn'].seqNum(): raise IOError('Failed to load sequences from seqin=FILE')
            return True     # Setup successful
        except: self.errorLog('Problem during %s setup.' % self.prog()); return False  # Setup failed
#########################################################################################################################
    def restSetup(self):    ### Sets up self.dict['Output'] and associated output options if appropriate.
        '''
        Run with &rest=docs for program documentation and options. A plain text version is accessed with &rest=help.
        &rest=OUTFMT can be used to retrieve individual parts of the output, matching the tabs in the default
        (&rest=format) output. Individual `OUTFMT` elements can also be parsed from the full (&rest=full) server output,
        which is formatted as follows:

        ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
        # OUTFMT:
        ... contents for OUTFMT section ...
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

        ### Available REST Outputs:
        There is currently no specific help available on REST output for this program.
        '''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for outfmt in self.restOutputOrder(): self.dict['Output'][outfmt] = 'No output generated.'
            #!# Add specific program output here. Point self.dict['Output'][&rest=X] to self.str key.
            return
        except: self.errorLog('RestSetup error')
#########################################################################################################################
    def restOutputOrder(self): return rje.sortKeys(self.dict['Output'])
#########################################################################################################################
    ### <3> ### Additional Class Methods                                                                                #
#########################################################################################################################
    def mtQuery(self): ### If circular, generate double-length mtDNA sequence
        '''
        If circular, generate double-length mtDNA sequence and update self.str['mtQuery']
        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not self.getBool('Circle'):
                self.setStr({'mtQuery':self.getStr('mtDNA')})
                self.printLog('#MTQRY','Using mtdna=FILE input {0} for mtDNA query (circle=F)'.format(self.getStr('mtDNA')))
                return True
            mt2x = rje.baseFile(self.getStr('mtDNA'),strip_path=True)+'2X.fasta'
            seq = self.obj['mtDNA'].seqs()[0]
            self.setStr({'mtQuery':mt2x})
            self.setInt({'mtLen':self.obj['mtDNA'].seqLen(seq)})
            self.printLog('#MTDNA','Mitochondrial DNA length: {0}'.format(rje_seqlist.dnaLen(self.getInt('mtLen'))))
            if rje.exists(mt2x) and not self.force():
                self.printLog('#MTQRY','Using existing {0} file for mtDNA query (force=F)'.format(mt2x))
                return True
            ### ~ [2] Generate double copy query ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            (seqname,sequence) = self.obj['mtDNA'].getSeq(seq)
            sname = self.obj['mtDNA'].shortName(seq)
            open(mt2x,'w').write('>{0}2X\n{1}{1}\n'.format(sname,sequence))
            self.printLog('#MTQRY','Output double sequence to {0} for mtDNA query (circle=T)'.format(mt2x))
            return True
        except: self.errorLog('%s.mtQuery error' % self.prog())
#########################################################################################################################
    def numtSearch(self):      ### Perform NUMT search (GABLAM by default)
        '''
        Perform NUMT search (GABLAM by default):
        +-- $BASEFILE.numtsearch.blast
        +-- $BASEFILE.numtsearch.gablam.tdt
        +-- $BASEFILE.numtsearch.hitsum.tdt
        +-- $BASEFILE.numtsearch.local.gff
        +-- $BASEFILE.numtsearch.local.tdt
        +-- $BASEFILE.numtsearch.log
        +-- $BASEFILE.numtsearch.unique.gff
        +-- $BASEFILE.numtsearch.unique.tdt
        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            wanted = ['.unique.tdt']
            gbase = '{0}.numtsearch'.format(self.basefile())
            gabrun = self.force() or not rje.checkForFiles(filelist=wanted,basename=gbase,log=self.log,cutshort=True,ioerror=False,missingtext='Not found.')
            if not gabrun:
                self.printLog('#SEARCH','NUMT Search results found (force=F). Will re-use')
                return True
            ### ~ [2] Perform NUMT search ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            gabdefaults = ['blastp=blastn','blasttask=blastn','blaste=1e-4','fasdir=numtfasta/','keepblast=T','fragfas=T']
            gabcmd = ['seqin={0}'.format(self.getStr('mtQuery')),'searchdb={0}'.format(self.getStr('SeqIn')),
                      'basefile={0}.numtsearch'.format(self.basefile()),'localgff=T','localsam=T',
                      'localunique=T','fullblast=T',
                      'localmin={0}'.format(self.getInt('MinFragLen'))]
            gabobj = gablam.GABLAM(self.log,gabdefaults+self.cmd_list+gabcmd)
            gabobj.run()
            rje.checkForFiles(filelist=wanted,basename=gbase,log=self.log,cutshort=True,ioerror=True,missingtext='Not found.')
            return
        except: self.errorLog('%s.numtSearch error' % self.prog())
#########################################################################################################################
    def numtProcess(self):      ### Load in NUMT results as numtfrag and merge to numtblock
        '''
        Load in NUMT results as numtfrag and merge to numtblock
        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            numtfrag = '{0}.numtsearch.unique.tdt'.format(self.basefile())
            mtlen = self.getInt('mtLen')
            ### ~ [2] Load NUMT fragments ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            fragdb = self.db().addTable(numtfrag,mainkeys=['Query','Hit','SbjStart','SbjEnd'],name='numtfrag',expect=True,replace=True,uselower=False)
            fragdb.dropFields(['Positives'])
            fragdb.dataFormat({'BitScore':'num','SbjStart':'int','SbjEnd':'int','QryStart':'int','QryEnd':'int','AlnID':'int','Expect':'num','Length':'int','Identity':'int'})
            if not fragdb: raise ValueError('Unable to load {0}!'.format(numtfrag))
            ## ~ [2a] Filter ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            filt = []
            mtid = max(0.0,self.getPerc('MTMaxID'))
            mtcov = max(0.0,self.getPerc('MTMaxCov'))
            for entry in fragdb.entries():
                if mtid and mtcov and ((entry['QryEnd'] - entry['QryStart'] + 1) / mtlen) > mtcov and (entry['Identity']/entry['Length']) > mtid:
                    filt.append(entry)
            if filt:
                filtseq = fragdb.dataList(filt,'Hit')
                fragdb.dropEntryList(filt,logtxt='Suspected non-NUMT mtDNA')
                fragdb.index('Hit')
                self.printLog('#MTFRAG','{0} suspected (non-NUMT) mtDNA sequence(s) and {1} fragment(s) (>{2:.1f}% coverage @ >{3:.1f}% identity)'.format(len(filtseq),len(filt),mtcov*100,mtid*100))
                for sname in filtseq:
                    if sname in self.list['Exclude']: self.printLog('#MTSEQ','{0} already in exclude=LIST'.format(sname))
                    elif self.getBool('MTMaxExclude'):
                        self.printLog('#MTSEQ','{0} added to exclude=LIST (mtmaxexclude=T)'.format(sname))
                        self.list['Exclude'].append(sname)
                    elif sname in fragdb.index('Hit'):
                        self.warnLog('Suspected mtDNA sequence {0} has additional unfiltered NUMT fragments'.format(sname))
            if self.list['Exclude']:
                fragdb.dropEntriesDirect('Hit',self.list['Exclude'])
            ## ~ [2b] Reformat ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            fragdb.renameField('Hit','SeqName')
            fragdb.renameField('SbjStart','Start')
            fragdb.renameField('SbjEnd','End')
            fragdb.renameField('QryStart','mtStart')
            fragdb.renameField('QryEnd','mtEnd')
            fragdb.addField('Strand',evalue='+')
            for entry in fragdb.entries():
                if entry['End'] < entry['Start']:
                    entry['Strand'] = '-'
                    (entry['End'],entry['Start']) = (entry['Start'],entry['End'])
            fragdb.setFields('SeqName	Start	End	Strand	BitScore	Expect	Length	Identity	mtStart	mtEnd'.split())
            fragdb.newKey(['SeqName','Start','End','Strand'])
            if self.getBool('Circle'):
                for entry in fragdb.entries():
                    if entry['mtStart'] > mtlen: entry['mtStart'] -= mtlen
                    if entry['mtEnd'] > mtlen: entry['mtEnd'] -= mtlen
            fragdb.saveToFile()

            ### ~ [3] Merge NUMT fragments ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            fragmerge = self.getInt('FragMerge')
            stranded = self.getBool('Stranded')
            blockdb = self.db().copyTable(fragdb,'numtblock')
            blockdb.addField('mtFrag',evalue='')
            blockdb.addField('FragNum',evalue=1)
            blockdb.addField('FragLen',evalue=0)
            blockdb.addField('FragGaps',evalue=0)
            prevfrag = None
            frags = blockdb.entries(sorted=True)
            blocks = []     # New list of entries for NUMT blocks
            while frags:
                self.progLog('\r#MERGE','Merging NUMT fragments within {0}: {1} frags -> {2} blocks.   '.format(rje_seqlist.dnaLen(fragmerge),rje.iLen(frags),rje.iLen(blocks)))
                nextfrag = frags.pop(0)
                nextfrag['FragLen'] = nextfrag['End'] - nextfrag['Start'] + 1
                merge = prevfrag and prevfrag['SeqName'] == nextfrag['SeqName'] and (nextfrag['Start'] - prevfrag['End']) <= fragmerge
                if stranded: merge = merge and prevfrag['Strand'] == nextfrag['Strand']
                if merge:
                    prevfrag['FragGaps'] += (nextfrag['Start'] - prevfrag['End'] - 1)
                    prevfrag['End'] = nextfrag['End']
                    prevfrag['Expect'] = min(nextfrag['Expect'],prevfrag['Expect'])
                    for field in ['BitScore','Length','Identity','FragNum','FragLen']:
                        prevfrag[field] += nextfrag[field]
                    prevfrag['mtFrag'] = '{0}|{1}-{2}'.format(prevfrag['mtFrag'],nextfrag['mtStart'], nextfrag['mtEnd'])
                    if prevfrag['Strand'] != nextfrag['Strand']: prevfrag['Strand'] = '+/-'
                else:
                    prevfrag = nextfrag
                    blocks.append(prevfrag)
                    prevfrag['mtFrag'] = '{0}-{1}'.format(prevfrag['mtStart'],prevfrag['mtEnd'])
                    continue
            self.printLog('\r#MERGE','Merging NUMT fragments within {0} complete: {1} frags -> {2} blocks'.format(rje_seqlist.dnaLen(fragmerge),rje.iStr(fragdb.entryNum()),rje.iLen(blocks)))
            ## ~ [3a] Update blockdb data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            blockdata = {}
            for entry in blocks:
                bkey = (entry['SeqName'],entry['Start'],entry['End'],entry['Strand'])
                entry['Length'] = entry['End'] - entry['Start'] + 1
                blockdata[bkey] = entry
            blockdb.dict['Data'] = blockdata
            blockdb.dropFields(['mtStart','mtEnd'])
            blockdb.saveToFile()
        except: self.errorLog('%s.numtProcess error' % self.prog()); raise
#########################################################################################################################
    def coverageOutputs(self):  ### Generates depth plots for mitochondrial genome
        '''
        Generates depth plots for mitochondrial genome.
        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            samcmd = ['seqin={0}'.format(self.getStr('mtQuery'))]
            sam = rje_samtools.SAMtools(self.log,['depthplot=T','readlen=T']+self.cmd_list+samcmd)
            sam.obj['DB'] = self.obj['DB']
            mtlen = self.getInt('mtLen')
            seq = self.obj['mtDNA'].seqs()[0]
            ### ~ [2] Generate 'rid' table from numtfrag table ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            fragdb = self.db('numtfrag')
            riddb = self.db().copyTable(fragdb, 'rid')
            #i# ['RID','Locus','Start','End','RLen','MLen','Clip5','Clip3']
            riddb.autoID()
            riddb.newKey(['#'])
            riddb.addField('Locus',evalue=self.obj['mtDNA'].shortName(seq))
            riddb.renameField('Length', 'RLen')
            riddb.setFields(['#','Locus','mtStart','mtEnd','RLen','Strand'])
            riddb.renameField('mtStart', 'Start')
            riddb.renameField('mtEnd', 'End')
            splitx = 0
            nextid = max(riddb.dataList(riddb.entries(),'#')) + 1
            self.printLog('#SPLIT','Splitting fragments that span circularisation for depth plots')
            for rentry in riddb.entries():
                #i# Split hits that go off the end
                if rentry['Strand'] == '+' and rentry['Start'] > rentry['End']:
                    newentry = rje.combineDict({},rentry)
                    newentry['End'] = mtlen
                    newentry['#'] = nextid
                    riddb.addEntry(newentry)
                    rentry['Start'] = 1
                    splitx += 1
                    nextid += 1
            self.printLog('#SPLIT','Split {0} NUMT fragments that span circularisation for depth plots'.format(splitx))
            riddb.renameField('#','RID')
            riddb.newKey(['RID'])
            if self.dev() or self.debugging(): riddb.saveToFile()
            ### ~ [3] Generate depth plot ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            sam.coverageFromRID(None, depthplot=sam.getBool('DepthPlot'))
            if sam.getBool('ReadLen'): sam.coverageFromRID(None, depthplot=True, readlen=True)
            return True
        except: self.errorLog('%s.coverageOutputs error' % self.prog()); return False
#########################################################################################################################
### End of SECTION II: New Class                                                                                        #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION III: MODULE METHODS                                                                                         #
#########################################################################################################################

#########################################################################################################################
### END OF SECTION III                                                                                                  #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION IV: MAIN PROGRAM                                                                                            #
#########################################################################################################################
def runMain():
    ### ~ [1] ~ Basic Setup of Program  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    try: (info,out,mainlog,cmd_list) = setupProgram()
    except SystemExit: return  
    except: rje.printf('Unexpected error during program setup:', sys.exc_info()[0]); return
    
    ### ~ [2] ~ Rest of Functionality... ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    try: NUMTFinder(mainlog,['basefile=numtfinder']+cmd_list).run()

    ### ~ [3] ~ End ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    except SystemExit: return  # Fork exit etc.
    except KeyboardInterrupt: mainlog.errorLog('User terminated.')
    except: mainlog.errorLog('Fatal error in main %s run.' % info.program)
    mainlog.endLog(info)
#########################################################################################################################
if __name__ == "__main__":      ### Call runMain 
    try: runMain()
    except: rje.printf('Cataclysmic run error: {0}'.format(sys.exc_info()[0]))
    sys.exit()
#########################################################################################################################
### END OF SECTION IV                                                                                                   #
#########################################################################################################################
