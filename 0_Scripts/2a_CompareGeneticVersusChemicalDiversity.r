#! /usr/bin/Rscript

# Compare hemp genetics and chemistry to see if more variable in one is more vairable in other

# Libraries
library(argparse)
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(vegan)

# Command-line arguments
options(stringsAsFactors=F)
parser=ArgumentParser()
parser$add_argument("-c", "--chemistry", help="Input file of HPLC results")
parser$add_argument("-d", "--distances", help="Input file of genetic distances (TASSEL format)")
parser$add_argument("-k", "--keyfile", help="Key file relating sample numbers to accessions")
parser$add_argument("-o", "--outprefix", help="Output file prefix")
parser$add_argument("-r", "--remove-outliers", help="Outliers to remove in second plot")
args=parser$parse_args()
# setwd('/home/jgwall/Projects/Hemp/HempAnalysis_JasonRedo/2_RevisionAdditions/')
# args=parser$parse_args(c("-c", "../0_HDGS_HPLC_Results.csv", "-d", "2a_hemp_distances.txt",
#                         "-k", "../0_Scripts/0_HDGS_keyfile.txt", "-o", "99_tmp", "-r", "Baox"))

cat("Plotting genetic versus chemical variance\n")

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
                                                     CBD_logvar = var(log(CBD+1)),
                                                     THC_logvar = var(log(THC+1)))

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
plot_graphic = function(mydata){
    label_key = c("CBD (raw)", "CBD (log)", "THC (raw)", "THC (log)")
    names(label_key) = c("CBD_var", "CBD_logvar", "THC_var", "THC_logvar")
    ggplot(mydata, mapping=aes(y=genetic_var, x=variance)) +
        geom_point() +
        geom_smooth(method=lm, fill='lightskyblue') +
        stat_cor(method = "pearson", label.y=-0.0002, size=3) +
        facet_wrap(~chemistry, scales="free_x", nrow=1, labeller=labeller(chemistry = label_key)) +
        labs(x="Phenotypic Variance", y="Average Genetic Distance")
}

allplot = plot_graphic(plotdata)

# Plot removing outliers
cat("\tRemoving outliers", args$remove_outliers, "\n")
cleandata =subset(plotdata, !plotdata$pop %in% args$remove_outliers) # Subset out the outliers
cleanplot = plot_graphic(cleandata)

# Write out
ggsave(allplot, file=paste(args$outprefix, "all.png", sep="."), width=7, height=3, units="in")
ggsave(cleanplot, file=paste(args$outprefix, "no_outliers.png", sep="."), width=7, height=3, units="in")
ggsave(allplot, file=paste(args$outprefix, "all.svg", sep="."), width=7, height=3, units="in")
ggsave(cleanplot, file=paste(args$outprefix, "no_outliers.svg", sep="."), width=7, height=3, units="in")
write.csv(variances, paste(args$outprefix, "all.csv", sep="."), row.names=FALSE)

################
# Check for individual plants
################

# Make distance matrices, setting 0 values to 1/10 of lowest non-zero value
chemistry = myphenos[,c('THC','CBD')]
chemistry$THC[chemistry$THC==0] = min(chemistry$THC[chemistry$THC>0]) * 0.1
chemistry$CBD[chemistry$CBD==0] = min(chemistry$CBD[chemistry$CBD>0]) * 0.1
chemistry = log10(chemistry) # Log of chemistry phenotypes
chemdist = dist(chemistry)
genetic.dist = as.dist(mydist)

# Confirm samples are in the same order
if(! identical(dimnames(as.matrix(chemdist)), dimnames(as.matrix(genetic.dist)))){
    stop("Chemical and genetic distance matrices do not match!!!")
}
mymantel = mantel(chemdist, genetic.dist, permutations=9999)
cat("Mantel test of genetic and chemical similarity has a p-value of", mymantel$signif, "\n")
