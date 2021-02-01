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
parser$add_argument("-o", "--outprefix", help="Output file prefix")
args=parser$parse_args()
# setwd('/home/jgwall/Projects/Hemp/CheckMattGenotyping/')
# args=parser$parse_args(c("-i", "1k_gs_results.txt", "-o", "99_tmp"))

# Load data
cat("Compiling and plotting genomic prediction accuracy\n")
preds = read.delim(args$infile)


# Calculate prediction accuracy
per_iteration = preds %>% group_by(Trait, Iteration) %>% summarize(iter_accuracy = mean(Accuracy))
stats = per_iteration %>% group_by(Trait) %>% summarize(mean_accuracy = mean(iter_accuracy), stdev_accuracy = sd(iter_accuracy))

# Graph
per_iteration$trait_type = sub(per_iteration$Trait, pattern="(THC|CBD|Ratio).+", repl="\\1")
myplot = ggplot(per_iteration, mapping=aes(x=Trait, y=iter_accuracy, fill=trait_type)) +
    geom_violin(draw_quantiles=c(0.5)) +
    theme(axis.text.x = element_text(angle=90)) +
    ggtitle("Genomic prediction accuracy")

# Write out results
write.table(stats, file=paste(args$outprefix, ".stats.txt",sep=''), sep='\t', quote=F, row.names=F, col.names=T)
ggsave(myplot, file=paste(args$outprefix, ".accuracy.png",sep=''))