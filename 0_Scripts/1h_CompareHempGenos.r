#! /usr/bin/Rscript

# Compare hemp trees and distances to see if close enough

# Libraries
library(argparse)
library(ggtree)
library(ggplot2)
library(gridExtra)
library(eulerr)
library(tidyverse)

# Command-line arguments
options(stringsAsFactors=F)
parser=ArgumentParser()
parser$add_argument("--oldtree", help="Newick-formatted tree for original genos")
parser$add_argument("--newtree", help="Newick-formatted tree for new genos")
parser$add_argument("--olddist", help="TASSEL-formatted distance matrix for original genos")
parser$add_argument("--newdist", help="TASSEL-formatted distance matrix for new genos")
parser$add_argument("--oldgenos", help="Hapmap-formatted original genos")
parser$add_argument("--newgenos", help="Hapmap-formatted new genos")
parser$add_argument("-k", "--keyfile", help="Keyfile of samples to population labels")
parser$add_argument("-o", "--outprefix", help="Output file prefix")
args=parser$parse_args()
setwd('/home/jgwall/Projects/Hemp/CheckMattGenotyping/')
# args=parser$parse_args(c("--oldtree", "1g_hemp_orig.tre",
#                          "--newtree", "1f_hemp_filtered.sorted.tre",
#                          "--olddist", "1g_hemp_orig.distances.txt",
#                          "--newdist", "1f_hemp_filtered.sorted.distances.txt",
#                          "--oldgenos", "../HempAnalysis_MattJohnson/HDGS_Analysis/Population_Structure/HDGS_less_het_no_Indel.hmp.txt",
#                          "--newgenos", "1f_hemp_filtered_.sorted.hmp.txt.gz",
#                          "-k", "0_HDGS_SplitsTree_KeyFile.txt",
#                          "-o", "99_tmp"))

# Load data
cat("Comparing old and new hemp genotypes\n")
oldtree = read.tree(args$oldtree)
newtree = read.tree(args$newtree)
olddist = read.delim(args$olddist, skip=5, header=F, row.names=1)
newdist = read.delim(args$newdist, skip=5, header=F, row.names=1)
oldgenos = read.delim(args$oldgenos, check.names=F)
newgenos = read.delim(args$newgenos, check.names=F)
key = read.delim(args$keyfile, row.names=1)
names(key)='pop'
het_calls = c("S","W", "K", "M", "Y","R")

# Reformat new labels to match old (New ones have some extra parts from pipeline)
rename_samples = function(x){
    sub(x, pattern=".+/(.+)\\.aligned.sorted.bam", repl="\\1")
}
rownames(newdist) = rename_samples(rownames(newdist))
colnames(newdist) = rownames(newdist)
colnames(olddist) = rownames(olddist)
newdist = newdist[rownames(olddist),colnames(olddist)] # Match order
newtree$tip.label = rename_samples(newtree$tip.label)
names(newgenos) = rename_samples(names(newgenos))

###################
# Raw genotypes
###################

# Check how many sites called in common
sites.old = paste(oldgenos$chrom, oldgenos$pos)
sites.new = paste(newgenos$chrom, newgenos$pos)
all_sites = unique(c(sites.old, sites.new))
sitekey = data.frame(row.names=all_sites, old_genos = all_sites %in% sites.old, new_genos = all_sites %in% sites.new)

png(paste(args$outprefix, ".sites.png", sep=""), width=600, height=600)
    plot(euler(sitekey), quantities = TRUE, main="Intersection of called sites")
dev.off()

# Check genotype concordance at common sites
common = intersect(sites.old, sites.new)
oldcalls = oldgenos[sites.old %in% common, 12:ncol(oldgenos)]
newcalls = newgenos[sites.new %in% common, 12:ncol(newgenos)]
newcalls = newcalls[,names(oldcalls)] # Match sample order
compare.calls = data.frame(old = unlist(oldcalls), new=unlist(newcalls)) %>%
    group_by(old, new) %>% summarize(count = length(old))

# Classify changes
compare.calls$type = "mismatch"
compare.calls$type[compare.calls$old == compare.calls$new] = "match"
compare.calls$type[(compare.calls$old %in% het_calls) | (compare.calls$new %in% het_calls)] = "change_het"
compare.calls$type[(compare.calls$old == "N") | (compare.calls$new=="N")] = "missing" # Be sure this comes last
compare.calls = compare.calls[order(compare.calls$type, compare.calls$count),]

# Plot genotype comparison
genoplot = ggplot(compare.calls, mapping=aes(x=type, y=count, fill=type)) +
    geom_bar(stat='sum') +
    theme(legend.position = 'none',
          axis.text.x = element_text(size=16)) +
    ggtitle(paste("Compare", length(unlist(oldcalls)), "genotypes at", length(common),"sites"))
ggsave(genoplot, file=paste(args$outprefix, ".genotypes.png", sep=""))

write.table(compare.calls, file=paste(args$outprefix, ".genotypes.txt", sep=""), quote=F, row.names=F, col.names=T, sep='\t')

###################
# Distance Matrices
###################

# Check that matrices match
if(!identical(rownames(newdist), rownames(olddist))){
    stop("Row names did not correctly match between new and old distance matrices")
}
if(!identical(colnames(newdist), colnames(olddist))){
    stop("Column names did not correctly match between new and old distance matrices")
}

# Assess correlation of distance matrices
newdist=as.dist(newdist)
olddist=as.dist(olddist)
mycor = cor(newdist, olddist)
dists = data.frame(old=as.numeric(olddist), new=as.numeric(newdist))
corplot = ggplot(dists, aes(x=old, y=new))+ 
    geom_point(alpha=0.1) +
    annotate("text", label=paste("r2 =", round(mycor, digits=3)), x=mean(olddist), y=max(newdist)) +
    ggtitle("Correlation of genetic distances")
ggsave(corplot, file=paste(args$outprefix, ".correlation.png", sep=""))

#################
# Phylogenetic Trees
#################


# Make tree keyfiles
make_treedata = function(mytree, key){
    treedata = ggtree(mytree)$data
    samplenum = sub(treedata$label, pattern="P.+-HDGS-", repl="")
    treedata$pop = key[samplenum, "pop"]
    return(treedata)
}
treedata.old = make_treedata(oldtree, key)
treedata.new = make_treedata(newtree, key)

# Helper function to plot mini trees below main one
minitree=function(mydata, color_column, colorkey){
    ggtree(mydata, mapping=aes(color=color_column), layout="fan", size=0.3) + 
        scale_colour_manual(values=colorkey)
}

# Make mini trees broken up by each population
make_minitrees=function(mypop, treedata, mycolor='purple'){
    is_target = treedata$pop == mypop
    minitree(treedata, color_column=is_target, colorkey=c("TRUE"=mycolor, 'FALSE'='gray70')) +
        ggtitle(mypop) +
        theme(plot.title = element_text(size=8, hjust=0.5, face='bold'))
}

pops = sort(unique(treedata.old$pop))
oldtree.plots = lapply(pops, make_minitrees, treedata=treedata.old, mycolor='red')
newtree.plots = lapply(pops, make_minitrees, treedata=treedata.new, mycolor='blue')

png(paste(args$outprefix, ".trees.png", sep=""), width=500, height=250*length(pops))
    grid.arrange(grobs=c(oldtree.plots, newtree.plots), ncol=2, as.table=FALSE)
dev.off()