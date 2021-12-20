## The Broad Institute
## SOFTWARE COPYRIGHT NOTICE AGREEMENT
## This software and its documentation are copyright (2012) by the
## Broad Institute/Massachusetts Institute of Technology. All rights are
## reserved.
##
## This software is supplied without any warranty or guaranteed support
## whatsoever. Neither the Broad Institute nor MIT can be responsible for its
## use, misuse, or functionality.

# Bridging script to accept parameters from GenePattern and use them to call HAPSEG

# Load required libraries
suppressMessages(suppressWarnings(library(boot)))
suppressMessages(suppressWarnings(library(class)))
suppressMessages(suppressWarnings(library(cluster)))
suppressMessages(suppressWarnings(library(foreign)))
suppressMessages(suppressWarnings(library(KernSmooth)))
suppressMessages(suppressWarnings(library(lattice)))
suppressMessages(suppressWarnings(library(MASS)))
suppressMessages(suppressWarnings(library(Matrix)))
suppressMessages(suppressWarnings(library(mgcv)))
suppressMessages(suppressWarnings(library(nlme)))
suppressMessages(suppressWarnings(library(nnet)))
suppressMessages(suppressWarnings(library(rpart)))
suppressMessages(suppressWarnings(library(spatial)))
suppressMessages(suppressWarnings(library(Rcpp)))
suppressMessages(suppressWarnings(library(numDeriv)))
suppressMessages(suppressWarnings(library(iterators)))
suppressMessages(suppressWarnings(library(foreach)))
suppressMessages(suppressWarnings(library(xtable)))
suppressMessages(suppressWarnings(library(getopt)))
suppressMessages(suppressWarnings(library(optparse)))
suppressMessages(suppressWarnings(library(DBI)))
suppressMessages(suppressWarnings(library(XML)))
suppressMessages(suppressWarnings(library(RSQLite)))
suppressMessages(suppressWarnings(library(RColorBrewer)))
suppressMessages(suppressWarnings(library(BiocGenerics)))
suppressMessages(suppressWarnings(library(IRanges)))
suppressMessages(suppressWarnings(library(Biobase)))
suppressMessages(suppressWarnings(library(AnnotationDbi)))
suppressMessages(suppressWarnings(library(annotate)))
suppressMessages(suppressWarnings(library(geneplotter)))
suppressMessages(suppressWarnings(library(DNAcopy)))
suppressMessages(suppressWarnings(library(HAPSEG)))
sessionInfo()

args <- commandArgs(trailingOnly=TRUE)

option_list <- list(
  make_option("--plate.name", dest="plate.name"),
  make_option("--array.name", dest="array.name"),
  make_option("--seg.file", dest="seg.file", default=NULL),
  make_option("--snp.file", dest="snp.file"),
  make_option("--genome.build", dest="genome.build"),
  make_option("--platform", dest="platform"),
  make_option("--use.pop", dest="use.pop"),
  make_option("--impute.gt", dest="impute.gt", type="logical"),
  make_option("--plot.segfit", dest="plot.segfit", type="logical"),
  make_option("--merge.small", dest="merge.small", type="logical"),
  make_option("--merge.close", dest="merge.close", type="logical"),
  make_option("--min.seg.size", dest="min.seg.size"),
  make_option("--is.normal", dest="normal", type="logical"),
  make_option("--out.p", dest="out.p"),
  make_option("--seg.merge.thresh", dest="seg.merge.thresh"),
  make_option("--use.normal", dest="use.normal", type="logical"),
  make_option("--adj.atten", dest="adj.atten", type="logical"),
  make_option("--drop.x", dest="drop.x", type="logical"),
  make_option("--drop.y", dest="drop.y", type="logical"),
  make_option("--calls.file", dest="calls.file", default=NULL),
  make_option("--mn.sample", dest="mn.sample", default=NULL),
  make_option("--out.file", dest="out.file", default=NULL),
  make_option("--calibrate.data", dest="calibrate.data", type="logical", default=NULL),
  make_option("--clusters.file", dest="clusters.file", default=NULL),
  make_option("--prev.theta.file", dest="prev.theta.file", default=NULL)
  )

opt <- parse_args(OptionParser(option_list=option_list), positional_arguments=TRUE, args=args)
print(opt)
opts <- opt$options

min.seg.size <- as.numeric(opts$min.seg.size)
if (is.na(min.seg.size)) {
  stop("min.seg.size must be a numeric value")
}

out.p <- as.numeric(opts$out.p)
if (is.na(out.p)) {
  stop("out.p must be a numeric value")
}

seg.merge.thresh <- as.numeric(opts$seg.merge.thresh)
if (is.na(seg.merge.thresh)) {
  stop("seg.merge.thresh must be a numeric value")
}

# returns string w/o leading or trailing whitespace
trim <- function (x) gsub("^\\s+|\\s+$", "", x)

use.normal <- opts$use.normal

## For production purposes, force.diploid should always be the same
## as normal.
force.diploid <- opts$normal

## If mn.sample is an empty string, set it to NULL
if (is.null(opts$mn.sample)) {
  mn.sample <- NULL
} else {
  mn.sample <- trim(opts$mn.sample)
  if (mn.sample == "") {
    if (use.normals) {
      stop("Please provide the matched normal sample name in mn.sample for use.normals=TRUE")
    }
    mn.sample <- NULL
  }
}

if (is.null(opts$calls.file)) {
  calls.file <- NULL
} else {
  calls.file <- trim(opts$calls.file)
  if (calls.file == "") {
    if (use.normals) {
      stop("A file is required for calls.file for use.normals=TRUE")
    }
    calls.file <- NULL
  } else {
    # Not yet sure whether this applies for standard end-users.  Probably not, but
    # we may have to account for the same thing somehow.  Could use a preprocessor
    # module ahead of this one or could add a flag to activate this code.
    #tmp.calls = read.delim(opt$calls.file)
    #tmp.calls = tmp.calls[-1, -3]
    #colnames(tmp.calls)[1] = "probeset_id"
    #calls.file = tempfile()
    #write.table(tmp.calls, file=calls.file, quote=FALSE, sep="\t", row.names=FALSE)
  }
}

if (is.na(opts$calibrate.data) || is.null(opts$calibrate.data)) {
  calibrate.data <- NULL
} else {
  calibrate.data <- opts$calibrate.data
}

if (is.null(opts$clusters.file)) {
  clusters.file <- NULL
} else {
  clusters.file <- trim(opts$clusters.file)
  if (clusters.file == "") {
    clusters.file <- NULL
  }
}

if (is.null(opts$prev.theta.file)) {
  prev.theta.file <- NULL
} else {
  prev.theta.file <- trim(opts$prev.theta.file)
  if (prev.theta.file == "") {
    prev.theta.file <- NULL
  }
}

results.dir <- getwd();

registerDoSEQ()

# The user will select a genome build via drop down, which will then pass a directory path 
# containing the relevant BEAGLE files.  The last piece of this path will be the name of the
# genome build. Only hg19 and hg18 are supported.
phased.bgl.dir <- opts$genome.build
genome.build <- basename(opts$genome.build)
if (!(genome.build %in% c("hg19", "hg18"))) {
   stop(paste0("Unrecognized genome build '", genome.build, "'"))
}

suppressWarnings(  # Done at JGentry's suggestion; R pkgs tend to put irrelevant output on stderr.

  RunHapSeg(opts$plate.name, opts$array.name, opts$seg.file,
            opts$snp.file, genome.build, results.dir,
            opts$platform, opts$use.pop, opts$impute.gt, opts$plot.segfit,
            opts$merge.small, opts$merge.close, min.seg.size,
            opts$normal, out.p, seg.merge.thresh,
            opts$use.normal, opts$adj.atten, phased.bgl.dir,
            force.diploid=force.diploid, drop.x=opts$drop.x, 
            drop.y=opts$drop.y, calls.fn=calls.file, 
            mn.sample=mn.sample, out.file=opts$out.file, 
            calibrate.data=calibrate.data, clusters.fn=clusters.file,
            prev.theta.fn=prev.theta.file, verbose=TRUE)
  )

sessionInfo()