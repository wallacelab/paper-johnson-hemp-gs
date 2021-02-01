#! /usr/bin/Rscript

# Check het calls to see if I trust the algorithm that called them

# Libraries
library(argparse)
library(ggplot2)
library(gridExtra)

# Command-line arguments
options(stringsAsFactors=F)
parser=ArgumentParser()
parser$add_argument("-i", "--infile", help="3-column input file of genotype calls, ref, and alt depths")
parser$add_argument("-p", "--prob-cutoff", type='double', default=0.01, help="3-column input file of genotype calls, ref, and alt depths")
parser$add_argument("-o", "--outfile", help="Output graphic file")
parser$add_argument("-l", "--log-scale", default=FALSE, action="store_true", help="Whether to log-scale the histogram counts")
args=parser$parse_args()
# setwd('/home/jgwall/Projects/Hemp/HempAnalysis_JasonRedo/1a_filtering/')
# args=parser$parse_args(c("-i","1l_check_het_calls.txt", "-o", "99_tmp.png"))


# Load data
cat("Graphing het call properties from",args$infile,"\n")
hets = read.delim(args$infile)

# Make sure just have het calls
is_het = hets$call == "0/1" | hets$call == "1/0"
if(any(!is_het)){
    cat("\tNon-het calls detected. Filtering.\n")
    hets = subset(hets, is_het)
    cat("\t\t",length(is_het),"total calls filtered down to", nrow(hets),"het calls\n")
}

# Cut down to just the allele depths
hets = hets[,c("ref","alt")]

# Get relative proportion of minor alleles
minor_count = apply(hets, MARGIN=1, FUN=min)
total_count = rowSums(hets)
minor_fraction = minor_count / total_count
ref_fraction = hets$ref / total_count
alt_fraction=  hets$alt / total_count

minorplot = qplot(x=minor_fraction, geom="histogram") +
    labs(x="Fraction reads", y="# Sites", title="Minor Allele") +
    geom_histogram(fill='cornflowerblue')

refplot = qplot(x=ref_fraction, geom="histogram") +
    labs(x="Fraction reads", y="# Sites", title="Reference Allele") +
    geom_histogram(fill='darkseagreen')

altplot = qplot(x=alt_fraction, geom="histogram") +
    labs(x="Fraction reads", y="# Sites", title="Alternate Allele") +
    geom_histogram(fill='khaki')

# Calculate p-values based on binomial distribution, assuming 50% probability for a het. 
# Basically, this is asking what's the p-value of getting this skewed a distribution _or more skewed_?
probs = pbinom(minor_count, size=total_count, prob=0.5) * 2    # * 2 because asking if is more skewed in either direction

# Note: For ones where it's a 50-50 split, you can actually get p-values >1 because the 
#   distributions "share" the 50/50 outcome, so it gets added twice. Turned these to just p=1
probs[probs>1] = 1

# Plot probabilities
too_small = sum(probs < args$prob_cutoff)/length(probs)
too_small = round(too_small*100, digits=1)
probplot = qplot(x=probs, geom='histogram')+
    labs(x="P-value from binomial", y="# Sites", title="P-values") +
    geom_histogram(fill='coral') +
    labs(subtitle=paste(too_small, "% are below ",args$prob_cutoff, sep=''))


# Assemble output
plots = list(minorplot, refplot, altplot, probplot)
if(args$log_scale){
    cat("Transforming count axis to Log10\n")
    for(i in 1:length(plots)){
        plots[[i]] = plots[[i]] + 
            scale_y_log10() +
            labs(y="log10(# sites)")
    }
}


png(args$outfile, width=600, height=1800)
    grid.arrange(grobs=plots, ncol=1)
dev.off()