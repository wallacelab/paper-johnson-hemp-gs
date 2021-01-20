#! /usr/bin/Rscript 

# Rename samples so can tell how close the clones are, etc.
library(argparse)

# Command-line arguments
options(stringsAsFactors=F)
parser=ArgumentParser()
parser$add_argument("-i", "--infile", help="Hapmap-formatted file of genotypes")
parser$add_argument("-k", "--keyfile", help="Keyfile of samples to population labels")
parser$add_argument("-o", "--outfile", help="Output file ")
args=parser$parse_args()
# setwd('/home/jgwall/Projects/Hemp/CheckMattGenotyping/')
# args=parser$parse_args(c("-i", "1f_hemp_filtered_.sorted.hmp.txt.gz", "-k", "0_HDGS_SplitsTree_KeyFile.txt","-o", "99_tmp.hmp.txt"))

# Load data
genos=read.delim(args$infile, check.names=F)
key=read.delim(args$keyfile, row.names=1)

# Correct genotype names to strip file data
names(genos) = sub(names(genos), pattern=".+HDGS-(.+)\\.aligned.sorted.bam", repl="\\1")

# Separate out metadata and actual calls
metadata=genos[,1:11]
calls = genos[,12:ncol(genos)]

# Match samples to pop names
popnames = key[names(calls),1]
names(calls) = paste(popnames, names(calls), sep="__")
names(calls) = gsub(names(calls), pattern=" ", repl="_")

# Rejoin and write out
renamed = cbind(metadata, calls)
write.table(renamed, file=args$outfile, sep='\t', quote=F, row.names=F, col.names=T)
