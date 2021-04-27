#! /usr/bin/Rscript

# Compare hemp genetics and chemistry to see if more vairable in one is more vairable in other

# Libraries
library(argparse)
library(ggplot2)
library(ggpubr)
library(tidyverse)

# Command-line arguments
options(stringsAsFactors=F)
parser=ArgumentParser()
parser$add_argument("-c", "--chemistry", help="Input file of HPLC results")
parser$add_argument("-d", "--distances", help="Input file of genetic distances (TASSEL format)")
parser$add_argument("-k", "--keyfile", help="Key file relating sample numbers to accessions")
parser$add_argument("-o", "--outprefix", help="Output file prefix")
args=parser$parse_args()
# setwd('/home/jgwall/Projects/Hemp/HempAnalysis_JasonRedo/2_RevisionAdditions/')
 # args=parser$parse_args(c("-c", "../0_HDGS_HPLC_Results.csv", "-d", "2a_hemp_distances.txt", 
 #                          "-k", "../0_Scripts/0_HDGS_keyfile.txt", "-o", "99_tmp"))

# Load data
phenos = read.csv(args$chemistry, row.names="Unique_ID")
distances = read.delim(args$distances, row.names=1, skip=5, header=F)
keyfile=read.delim(args$keyfile)

# Rename samples so completely match
rownames(phenos) = toupper(rownames(phenos))
rownames(distances) = sub(rownames(distances), pattern="1_bams_sorted/P.-HDGS-(.+).aligned.sorted.bam", repl="HDGS_\\1")
colnames(distances) = rownames(distances)
key = keyfile$Populations
names(key) = paste("HDGS_", keyfile$sample, sep="")


# Subset so completely match
taxa = intersect(rownames(phenos), rownames(distances))
mydist = distances[taxa, taxa]
myphenos = phenos[taxa,]
mykey = key[taxa]

# Get chemical variance by accession
myphenos$pop = key[rownames(myphenos)]
variances = myphenos %>% group_by(pop) %>% summarize(CBD_var=var(CBD), THC_var=var(THC),
                                                     CBD_logvar = var(log(CBD)),
                                                     THC_logvar = var(log(THC)))

# Get genetic distances within each accession
genetic_var = sapply(unique(mykey), function(accession){
    targets = names(mykey)[mykey==accession]
    target_dists = mydist[targets, targets]  # Subset out distance matrix
    target_dists = as.numeric(as.dist(target_dists)) # Convert to non-redundant vector vector
    var(target_dists) # Return variances in genetic distances
})


# Combine genetic variance with chemical
variances$genetic_var = genetic_var[variances$pop]

# Plot
plotdata = pivot_longer(variances, cols=c("CBD_var", "CBD_logvar", "THC_var", "THC_logvar"), 
                        names_to="chemistry", values_to="variance")
ggplot(plotdata, mapping=aes(x=genetic_var, y=variance)) +
    geom_point() +
    geom_smooth(method=lm) +
    stat_cor(method = "pearson") +
    facet_wrap(~chemistry, scales="free_y")

# Write out
ggsave(paste(args$outprefix, ".png", sep=""))
write.csv(variances, paste(args$outprefix, ".csv", sep=""), row.names=FALSE)
