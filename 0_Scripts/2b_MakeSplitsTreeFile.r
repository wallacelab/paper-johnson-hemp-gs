#! /usr/bin/Rscript

# Compare hemp genetics and chemistry to see if more variable in one is more vairable in other

# Libraries
library(argparse)
library(RSplitsTree)

# Command-line arguments
options(stringsAsFactors=F)
parser=ArgumentParser()
parser$add_argument("-d", "--distances", help="Input file of genetic distances (TASSEL format)")
parser$add_argument("-k", "--keyfile", help="Key file relating sample numbers to accessions")
parser$add_argument("-o", "--outfile", help="Output file (nexus format)")
args=parser$parse_args()
# setwd('/home/jgwall/Projects/Hemp/HempAnalysis_JasonRedo/2_RevisionAdditions/')
# args=parser$parse_args(c("-d", "2a_hemp_distances.txt",
#                         "-k", "../0_Scripts/0_HDGS_keyfile.txt", "-o", "99_tmp.nex"))

cat("Converting distances to Nexus format\n")

# Load data
distances = read.delim(args$distances, row.names=1, skip=5, header=F)
keyfile=read.delim(args$keyfile)

# Format distance matrix
rownames(distances) = sub(rownames(distances), pattern="1_bams_sorted/P.-HDGS-(.+).aligned.sorted.bam", repl="HDGS_\\1")
distances = as.matrix(distances)

# Replace labels with names
rownames(keyfile) = paste("HDGS_", keyfile$sample, sep="")
keyfile$Populations = paste(keyfile$Populations, keyfile$sample, sep=".")
rownames(distances) = keyfile[rownames(distances),"Populations"]
distances = as.dist(distances)

# Convert to nexus
splitstree(distances, nexus.file=args$outfile)
