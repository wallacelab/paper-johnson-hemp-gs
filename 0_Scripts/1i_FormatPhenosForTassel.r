#! /usr/bin/Rscript

# Compare hemp trees and distances to see if close enough

# Libraries
library(argparse)
library(ggplot2)
library(tidyverse)

# Command-line arguments
options(stringsAsFactors=F)
parser=ArgumentParser()
parser$add_argument("-i", "--infile", help="Input file of phenotypes")
parser$add_argument("-o", "--outprefix", help="Output file prefix for graphic and TASSEL-formatted phenotypes")
parser$add_argument("-l", "--log", default=FALSE, action="store_true", help="Whether to add log+1-transformed phenotypes")
parser$add_argument("-s", "--sqrt", default=FALSE, action="store_true", help="Whether to add square root-transformed phenotypes")
args=parser$parse_args()
# setwd('/home/jgwall/Projects/Hemp/CheckMattGenotyping/')
# args=parser$parse_args(c("-i", "0_HDGS_HPLC_Results.csv", "-o", "99_tmp", "--log", "--sqrt"))

# Load data
cat("Transforming phenotypes. WARNING: Phenotype names, etc. are hardcoded, so be careful when adapting")
phenos = read.csv(args$infile, row.names=1, na.strings=c("NA", "na"))

# Log-transform
if(args$log){
    cat("\tLog-transforming CBD and THC, and getting a new ratio of their transformed values\n")
    phenos$THC_log = log(phenos$THC+1)
    phenos$CBD_log = log(phenos$CBD+1)
    phenos$Ratio_of_logs = phenos$CBD_log / phenos$THC_log
    phenos$Ratio_of_logs[is.infinite(phenos$Ratio_of_logs)] = NA
}

# Square root-transform
if(args$sqrt){
    cat("\tSquare root-transforming CBD and THC, and getting a new ratio of their transformed values\n")
    phenos$THC_sqrt = sqrt(phenos$THC)
    phenos$CBD_sqrt = sqrt(phenos$CBD)
    phenos$Ratio_of_sqrt = phenos$CBD_sqrt / phenos$THC_sqrt
    phenos$Ratio_of_sqrt[is.infinite(phenos$Ratio_of_sqrt)] = NA
}

# Format tidily for graphing
plotdata = pivot_longer(phenos, cols = names(phenos), names_to = "pheno")
plotdata$trait = sub(plotdata$pheno, pattern="(THC|CBD|Ratio).+", repl="\\1")
plotdata$trans = "none"
plotdata$trans[grep(plotdata$pheno, pattern="log")] = "log"
plotdata$trans[grep(plotdata$pheno, pattern="sqrt")] = "sqrt"

# Graph
plots = ggplot(plotdata, mapping=aes(x=value, fill=trait)) +
    geom_histogram() +
    facet_wrap(trait ~ trans, scales="free")

# Output Graphic
ggsave(plots, file=paste(args$outprefix, ".distributions.png", sep=""))

# Write out
outfile = paste(args$outprefix, ".tassel.txt", sep="")
write("<Phenotype>", file=outfile)
write(c("taxa", rep("data", ncol(phenos))), file=outfile, append=T, ncol=ncol(phenos)+1, sep='\t')
output = data.frame(rownames(phenos), phenos)
names(output)[1]="Taxon"
write.table(output, file=outfile, sep='\t', quote=F, row.names=F, col.names=T, append=T)



