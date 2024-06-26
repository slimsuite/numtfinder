########################################################
### SAMTools analysis data visualisation script      ###
### VERSION: 0.3.0                             ~~~~~ ###
### LAST EDIT: 16/06/17                        ~~~~~ ###
### AUTHORS: Richard Edwards 2016              ~~~~~ ###
### CONTACT: richard.edwards@unsw.edu.au       ~~~~~ ###
########################################################

# This script is for plotting sequencing coverage depth. The main functions for doing so will be
# placed in rje_genomics.r so that they are available for other scripts, e.g. PAGSAT.R.

# The required input files for these "samdepth" plots are (where * is basefile):
# *.depthplot.tdt
# *.coverage.tdt
# Optionally: *.check.tdt

# See notes below for "samunzip" SAMPhaser plots.


################# ::: HISTORY ::: ######################
# v0.0.0 : Initial version based on PAGSAT v0.7.4.
# v0.1.0 : Added SAMPhaser methods for used with SAMPhaser v0.4.0.
# v0.2.0 : Added read length plots.
# v0.3.0 : Added dirlen plot data where available.


#!# Rearrange a little and make more elements into functions (load data etc.)

############### ::: GENERAL SETUP ::: ##################
# Usage = Rscript rje.r samtools basefile [options] 

setTesting = function(){
  arglist = list(rdir=rdir,rtype="samdepth",basefile="/Users/redwards/Documents/Projects/Phasing-May16/dev/2016-06-21.S288C.TestSAM/S288C-ISH.chrI.depthtest")
  #arglist$basefile = "/Users/redwards/Data/projects/YeastPacBio-Jun15/analysis/2015-08-18.PAGSAT.MBG479/MBG479.SP16508.hcq.qv20.sgd.srt.PAGSAT/MBG479.SP16508.hcq.qv20.sgd.srt.L500ID800"
}
#i# Run setTesting() if testing within RStudio for development.

## ~ Read commandline arguments ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
rtype = arglist$rtype        # arglist is generated by rje.r or setTesting()
basefile = arglist$basefile  # arglist is generated by rje.r or setTesting()
if(rtype == "samunzip"){
  sambase = arglist$sambase
}else{
  sambase = basefile
}
print(basefile)

## ~ Set up global parameters and file names ~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
settings = list(basefile=basefile)   # Access with settings[["basefile"]] or settings$basefile
# General output options
settings$pngwidth = 2400
settings$pngheight = 1600
settings$pointsize = 36
# Genome plot options
settings$units = "kb"
# SAMTools depth plot options 
settings$logplot = TRUE  # Whether to plot Xdepth on log scale

# OLD PAGSAT options not used >>>
#settings$refbase = refbase 
#settings$assessment = TRUE # Process outputs from assessment run
#settings$covplot = TRUE    # Whether to generate chromosome and contig coverage plots.
#settings$diploid = TRUE    # Whether to run ChromAlignLoc in diploid mode.
#settings$genesummary = TRUE    # Whether to run GeneSummary plot.
#settings$chromalign = TRUE    # Whether to run ChromAlign plot.
#settings$protsummary = TRUE    # Whether to run ProtSummary plot.
#settings$features = TRUE    # Whether to expect Features table plot.
#settings$minloclen = 1000
#settings$topcontigs = 8        # Number of contigs to display as local alignment matches
#settings$topcontigprop = 0.33  # Min proportion of contig that must match chromosome for local alignment
settings$plotft = c("gene","mRNA","CDS","rRNA","tRNA","ncRNA","mobile","LTR","origin","centromere","telomere")
# <<< Delete when save to do so.

## ~ Update default settings from commandline arguments (setting=value) ~~~~~ ##
# List of boolean parameters:
boollist = c("logplot")
# List of numerical parameters:
numlist = c("pngwidth","pngheight","pointsize")
# List of list parameters:
listlist = c()
# Update arglist and then settings:
(arglist)  # This contains settings read in from arguments
arglist = argBool(arglist,boollist)    # List of Boolean settings
arglist = argNum(arglist,numlist)     # List of numeric settings
arglist = argList(arglist,listlist)
for(cmd in names(arglist)){
  settings[[cmd]] = arglist[[cmd]]  # Need to add processing of booleans etc.
}
(settings)  # This contains settings read in from arguments

## ~ Setup base filenames for outputs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
settings$basename = strsplit(basefile,'/')[[1]]
(settings$basename = settings$basename[length(settings$basename)])
plotdir = "SAMPlots"
if(rtype=="samunzip"){ plotdir = "haplotigs.Plots" }
plotdir = paste(basefile,".",plotdir,"/",sep="")
(settings$pngbase = paste(plotdir,settings$basename,sep=""))
if(file.exists(plotdir)==FALSE){
  dir.create(plotdir)
}

## ~ Setup scaling and units ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
#?# Make an rje_genomics function?
settings$scaling = 1
if(settings$units == "kb"){
  settings$scaling = 1000
}
if(settings$units == "Mb"){
  settings$scaling = 1e6
}

## ~ Setup colour profile ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
settings$col = list()   # This will have contigs and chromosomes added later.
# Setup feature colour scheme
#for(fi in 1:length(settings$plotft)){
#  settings$col[[settings$plotft[fi]]] = soton$col[fi+1]
#}
#settings$col[["LTR"]] = soton$col[15]
#settings$col[["mobile"]] = soton$col[16]

(settings)  # This contains final settings


############# ::: LOAD GENOME PLOTTING FUNCTIONS ::: ##################
rjesource("rje_genomics.r")


############# ::: LOAD AND TIDY DATA ::: ##################
#i# This script is now used for multiple SAM plotting functions. They will not all have the same input files.

#># rtype == "samdepth"
  # |-- *.depthplot.tdt : depthdb$Basefile = Main plot data
  # |-- *.coverage.tdt  : xcovdb$Basefile = Used for median XCoverage and Locus names
  # |-- [*.check.tdt]
#># rtype == "samreadlen" -> Like the depth plot but with max read length instead of Xdepth
# |-- *.readlenplot.tdt : depthdb$Basefile = Main plot data
# |-- *.readlen.tdt  : xcovdb$Basefile = Used for median XCoverage and Locus names
# |-- *.dirnlenplot.tdt : depthdb$DirnLen = Directional read length data
#># rtype == "samunzip" -> Stack plot of each starting locus and haplotigs produced. [basefile.locus.haplotigs.png]
  # |-- *.depthplot.tdt : depthdb$Basefile
  # |-- *.coverage.tdt  : xcovdb$Basefile
  # |-- *.haplotigs.depthplot.tdt 
  # |-- *.haplotigs.coverage.tdt
  # |-- *.haplotigs.tdt = positions of haplotigs in main loci

## ~ Detailed depthplot data ~ ##
depfile = list(samdepth="depthplot.tdt",samunzip="depthplot.tdt",samreadlen="readlenplot.tdt",sampart="depthplot.tdt")
depthdb = list()
# Three fields: Locus, Pos, X
dtdtfile = paste(sambase,depfile[[rtype]],sep=".")
if(! file.exists(dtdtfile)){
  dtdtfile = paste(basefile,depfile[[rtype]],sep=".")
}  
depthdb$Basefile = read.table( dtdtfile, header = T, stringsAsFactors = T, sep = '\t', quote = '', comment.char="")
#colnames(depthdb$Basefile) = c("Locus",'Pos','X')
summary(depthdb$Basefile)

## ~ Directional read length ~ ##
settings$dirndata = FALSE
if(rtype == "samreadlen"){
  xtype = "DirnLen"
  dfile = paste(sambase,"dirnlenplot.tdt",sep=".")
  if(file.exists(dfile)){
    depthdb[[xtype]] = read.table( dfile, header = T, stringsAsFactors = T, sep = '\t', quote = '', comment.char="")
    colnames(depthdb[[xtype]]) = c("Locus",'Pos','Len5','Len3')
    depthdb[[xtype]]$Len5 = depthdb[[xtype]]$Len5/1000
    depthdb[[xtype]]$Len3 = depthdb[[xtype]]$Len3/1000
    summary(depthdb[[xtype]])
    settings$dirndata = TRUE
  }
  
}


for(xtype in c("haplotigs")){
  dfile = paste(basefile,xtype,"depthplot.tdt",sep=".")
  if(file.exists(dfile)){
    depthdb[[xtype]] = read.table( dfile, header = T, stringsAsFactors = T, sep = '\t', quote = '', comment.char="")
    colnames(depthdb[[xtype]]) = c("Locus",'Pos','X')
    depthdb[[xtype]]$HapAcc = as.character(depthdb[[xtype]]$Locus)
    #!# Make this more efficient by cycling through Locus levels!
    for(i in 1:nrow(depthdb[[xtype]])){
      depthdb[[xtype]]$HapAcc[i] = strsplit(depthdb[[xtype]]$HapAcc[i],"__")[[1]][2]
    }
    summary(depthdb[[xtype]])
  }
  # Read readlenplot data for haplotigs if it exists
  dfile = paste(basefile,xtype,"readlenplot.tdt",sep=".")
  xtyper = paste(xtype,"readlen",sep=".")
  if(file.exists(dfile)){
    depthdb[[xtyper]] = read.table( dfile, header = T, stringsAsFactors = T, sep = '\t', quote = '', comment.char="")
    colnames(depthdb[[xtyper]]) = c("Locus",'Pos','X')
    summary(depthdb[[xtyper]])
  }
}

## ~ ReadLen Summary for SAMUnzip ~ ##
dfile = paste(sambase,depfile[["samreadlen"]],sep=".")
settings$readlenplot = FALSE
if(file.exists(dfile)){
  settings$readlenplot = TRUE
  depthdb$ReadLen = read.table( paste(sambase,depfile[["samreadlen"]],sep="."), header = T, stringsAsFactors = T, sep = '\t', quote = '', comment.char="")
  colnames(depthdb$ReadLen) = c("Locus",'Pos','X')
}

## ~ Coverage Summary ~ ##
xfile = list(samdepth="coverage.tdt",samunzip="coverage.tdt",samreadlen="readlen.tdt",sampart="coverage.tdt")
xcovdb = list()
# Should have MedianX
xcovfile = paste(sambase,xfile[[rtype]],sep=".")
if(! file.exists(xcovfile)){
  xcovfile = paste(basefile,xfile[[rtype]],sep=".")
}  
xcovdb$Basefile = read.table( xcovfile, header = T, stringsAsFactors = T, sep = '\t', quote = '', comment.char="")
summary(xcovdb$Basefile)
(medianx = median(xcovdb$Basefile$MedianX))

## ~ ReadLen Summary for SAMUnzip ~ ##
xcovfile = paste(sambase,xfile[["samreadlen"]],sep=".")
if(file.exists(xcovfile)){
  xcovdb$ReadLen = read.table( xcovfile, header = T, stringsAsFactors = T, sep = '\t', quote = '', comment.char="")
  (medianxr = median(xcovdb$ReadLen$MedianX))
}

## ~ Coverage Check Regions ~ ##
settings$checkfile = file.exists(paste(basefile,".check.tdt",sep=""))
if(rtype == "samreadlen"){ settings$checkfile==FALSE }
if(settings$checkfile==TRUE){
  checkdata = read.table( paste(basefile,".check.tdt",sep=""), header = T, stringsAsFactors = T, sep = '\t', quote = '', comment.char="")
  summary(checkdata)
  checkdesc = FALSE  # Whether description field included in file
  spans = c()
  for(field in colnames(checkdata)){
    if(substr(field,1,4) == "Desc"){
      checkdesc = TRUE
    }
    if(field == "Chrom"){
      checkdata$Locus = checkdata$Chrom
    }
    if(field == "Contig"){
      checkdata$Locus = checkdata$Contig
    }
    if(substr(field,1,4) == "Span"){
      spanx = as.integer(substr(field,5,nchar(field)))
      spans = c(spans,spanx)
    }
  }
}


### Haplotig Data
#i# hapdb$hapsnp = *.hapsnp.tdt
#># Locus	Block	Pos	Ref	A	B	pA
#i# hapdb$haplotigs = *.haplotigs.tdt
#># Locus	Block	Track	Start	End	HapStart	HapEnd	Haplotig	nSNP	nRID	MeanX
hapdb = list()
for(xtype in c("haplotigs","hapsnp")){
  dfile = paste(basefile,xtype,"tdt",sep=".")
  if(file.exists(dfile)){
    hapdb[[xtype]] = read.table( dfile, header = T, stringsAsFactors = T, sep = '\t', quote = '', comment.char="")
    summary(hapdb[[xtype]])
  }else{ hapdb[[xtype]] = data.frame()}
}

# Update colour dicionary for chromosomes (not used!)
cx = 0
for(chrom in xcovdb$Basefile[order(xcovdb$Basefile$Length,decreasing=TRUE),"Locus"]){
  cx = cx + 1
  if(cx > length(soton$col)){ cx = 1 }
  settings$col[[chrom]] = soton$col[cx]
  cname = strsplit(chrom,"_")[[1]][1]
  settings$col[[cname]] = soton$col[cx]
}


############# ::: BASIC SAMTOOLS DEPTHPLOT ::: ##################
if(rtype %in% c("samdepth","sampart")){
  # Generate plot for each chromosome
  for(chrom in levels(depthdb$Basefile$Locus)){
    pngfile=paste(settings$pngbase,"depth",chrom,"png",sep=".")
    png(filename=pngfile, width=settings$pngwidth, height=settings$pngheight, units = "px", pointsize=settings$pointsize)
    chromDepthPlot(depthdb$Basefile,chrom,medianx,partplot=rtype=="sampart")
    #!# Add checkdata to chromDepthPlot
    if(settings$checkfile==TRUE){
      regionDepth(checkdata,chrom,spans)
    }
    dev.off()
  }
}

############# ::: BASIC READLEN DEPTHPLOT ::: ##################
if(rtype=="samreadlen"){
  # Generate plot for each chromosome
  for(chrom in levels(depthdb$Basefile$Locus)){
    pngfile=paste(settings$pngbase,"readlen",chrom,"png",sep=".")
    png(filename=pngfile, width=settings$pngwidth, height=settings$pngheight, units = "px", pointsize=settings$pointsize)
    if(settings$dirndata == TRUE){
      chromDepthPlot(depthdb$Basefile,chrom,medianx,dirplot=TRUE,dirdata=depthdb$DirnLen,ylab="Read Len (kb)")
    }else{
      chromDepthPlot(depthdb$Basefile,chrom,medianx,ylab="Read Len (kb)")
    }
    #!# Add checkdata to chromDepthPlot
    if(settings$checkfile==TRUE){
      regionDepth(checkdata,chrom,spans)
    }
    dev.off()
  }
}


############# ::: SAMPHASER STACKED DEPTHPLOTS ::: ##################
#!# Can we make this a three panel plot with readlen data as well? #!#
#!# if(settings$readlenplot == TRUE){ }
# This will be a two panel plot:
# 1. Depth plot for original locus with SNP frequency and phased regions marked
#=# chromDepthPlot(chromdata,chrom,onex)
#=# regionDepth(checkdata,chrom) : Needs Locus, Start, End, Span0, MeanX
# 2. A triple stack of A, B and C haplotig regions based on the PAGSAT hitstack
#=# samphaserStackPlot(depdata,hapdepdata,hapdata,locus,onex)
#i# hapdb$hapsnp = *.hapsnp.tdt
#># Locus	Block	Pos	Ref	A	B	pA
#i# hapdb$haplotigs = *.haplotigs.tdt
#># Locus	Block	Track	Start	End	HapStart	HapEnd	Haplotig	nSNP	nRID	MeanX
if(rtype=="samunzip"){
  # Cycle through each locus
  for(locus in levels(hapdb$haplotigs$Locus)){
    # Setup PNG file
    pngfile=paste(settings$pngbase,"haplotigs",locus,"png",sep=".")
    png(filename=pngfile, width=settings$pngwidth, height=settings$pngheight, units = "px", pointsize=settings$pointsize)
    # Setup multipanel plot space
    #panels = c(rep(1,1),rep(2,3))  # Entire bottom panel to have room for contigs
    if(settings$readlenplot){
      layout(matrix(1:2,byrow=TRUE,nrow=2))
    }
    #par(mar=c(2,6,2,1))
    # 1. Depth plot for original locus with SNP frequency and phased regions marked
    #chromDepthPlot(depthdb$Basefile,locus,medianx)
    #!# This is now in samphaserStackPlot2()
    #i# regionDepth(checkdata,chrom) : Needs Locus, Start, End, Span0, MeanX
    # 2. A triple stack of A, B and C haplotig regions based on the PAGSAT hitstack
    par(mar=c(5,6,2,1))
    samphaserStackPlot2(depthdb$Basefile,depthdb$haplotigs,hapdb$haplotigs,hapdb$hapsnp,locus,medianx)
    if(settings$readlenplot){
      samphaserStackPlot2(depthdb$ReadLen,depthdb$haplotigs.readlen,hapdb$haplotigs,hapdb$hapsnp,locus,medianxr)
    }
    # Finish    
    dev.off()
  }
}

