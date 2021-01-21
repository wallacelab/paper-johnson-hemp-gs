#! /usr/bin/Rscript

# Test for distortion in the kmers made from each sequence

# Libraries
library(argparse)
library(ggplot2)
library(gridExtra)

# Command-line arguments
options(stringsAsFactors=F)
parser=ArgumentParser()
parser$add_argument("-s", "--sitefile", help="TASSEL genotype summary file of sites")
parser$add_argument("-t", "--taxafile", help="TASSEL genotype summary file of taxa")
parser$add_argument("-d", "--depthfile", help="VCFTools depth output file")
parser$add_argument("-g", "--genodepth", help="VCFTools output of genotype depth")
parser$add_argument("-r", "--readdepth", help="LGC genomics file of sample read depths")
parser$add_argument("-o", "--outfile", help="Output graphic")
parser$add_argument("-m", "--max-datapoints", type="integer", help="Randomly subset to this many datapoints if dataset has more")
args=parser$parse_args()
# setwd('/home/jgwall/Projects/Hemp/HempAnalysis_JasonRedo/1a_filtering/')
# args=parser$parse_args(c("-s","1c_sitesummary.txt", "-t", "1c_taxasummary.txt", 
#                          "-d", "1c_site_depth.txt", "-g", "1c_geno_depth.txt.gz", 
#                          "-r", "../0_sample_read_counts.tsv", "-o", "99_tmp.png"))


plots = list()


###
# Depth summary 
###


if(!is.null(args$depthfile)){
    
    cat("Processing site depth data from",args$depthfile, "\n")
    depthvals = read.delim(args$depthfile, check.names=F)$SUM_DEPTH #read in 
    if(!is.null(args$max_datapoints) && (length(depthvals) > args$max_datapoints)){ # Subset if needed
        cat("\tRandomly subsetting to",args$max_datapoints,"data points\n")
        depthvals = sample(depthvals, args$max_datapoints, replace=F)
    }
    depthvals = sort(depthvals) # Sort for plotting
    
    # Depth - TODO: As histogram or as dot plot? (Prob dots to match % het)
    depth_raw = qplot(1:length(depthvals), depthvals, color=I('darkgoldenrod')) +
        labs(title = "Sites - Read Depth", x="Site # (sorted)", y="Read Depth") 

    depth_log = qplot(1:length(depthvals), log10(depthvals), color=I('darkgoldenrod4')) +
        labs(title = "Sites - Read Depth (Log scale)", x="Site # (sorted)", y="Log10 Read Depth") 
    
    
    plots = c(plots, list(depth_raw, depth_log))
    
}

if(!is.null(args$genodepth)){
    
    cat("Processing genotype depth data from",args$genodepth, "\n")
    cat("\t(Note: this could take a while)\n")
    
    # Use header to determine number of samples to speed reading in
    header=scan(args$genodepth, nline=1, what=character())
    ntaxa = length(header)-2 # Take off 2 for CHROM, POS 
    
    # Read in file
    colClasses = c("NULL", "NULL", rep("numeric", ntaxa)) #Specifying column classes speeds things up
    genodepth = read.delim(args$genodepth, header=TRUE, colClasses = colClasses)
    genodepth = as.matrix(genodepth) # Also to speed things up
    colnames(genodepth) = sub(colnames(genodepth), pattern=".+(P..HDGS.[0-9]+).+", repl="\\1")
    colnames(genodepth) = gsub(colnames(genodepth), pattern=".", repl="-", fixed=T)
    
    # Normalize to read depth per sample
    sample_depths = read.delim(args$readdepth, row.names=1)
    read_depth = sample_depths[colnames(genodepth),"Restriction.enzyme.filtered.read.pairs"]
    norm_factors = mean(read_depth) / read_depth
    genodepth_norm = t(t(genodepth) * norm_factors)
    
    # Binplot of depths
    depthstats = data.frame(total_depth = rowSums(genodepth),
                            mean_depth = apply(genodepth, FUN=mean, MARGIN=1),
                            sd_depth = apply(genodepth, FUN=sd, MARGIN=1),
                            total_depth_norm = rowSums(genodepth_norm),
                            mean_depth_norm = apply(genodepth_norm, FUN=mean, MARGIN=1),
                            sd_depth_norm = apply(genodepth_norm, FUN=sd, MARGIN=1),
                            n_present = rowSums(genodepth > 0))
    depthstats$fraction_present = depthstats$n_present / ntaxa
    
    ## Binplot - Everything
    bindepth = ggplot(depthstats, mapping=aes(x=log10(mean_depth_norm), y=fraction_present)) +
        geom_bin2d() +
        ggtitle("Mean depth vs faction present (all)") +
        scale_fill_gradient(low='blue', high='red')
    
    # Binplot - Only things with at least 0.5 mean depth
    substats = subset(depthstats, depthstats$mean_depth >= 0.5)
    ggplot(substats, mapping=aes(x=log10(mean_depth_norm), y=fraction_present)) +
        geom_bin2d() +
        ggtitle("Mean depth vs faction present (mean >= 0.5)")+
        scale_fill_gradient(low='blue', high='red')
    
    
    if(!is.null(args$max_datapoints) && length(genodepth) > args$max_datapoints){ # Subset if needed
        cat("\tRandomly subsetting to",args$max_datapoints,"data points\n")
        genodepth = sample(genodepth, args$max_datapoints, replace=F)
    }
    genodepth = sort(genodepth) # Sort for plotting
    
    
    # Depth - TODO: As histogram or as dot plot? (Prob dots to match % het)
    genodepth_raw = qplot(1:length(genodepth), genodepth, color=I('deeppink')) +
        labs(title = "Geno Calls - Read Depth", x="Site # (sorted)", y="Read Depth") 

    genodepth_log = qplot(1:length(genodepth), log10(genodepth), color=I('deeppink4')) +
        labs(title = "Geno Calls - Read Depth (Log scale)", x="Site # (sorted)", y="Log10 Read Depth") 
    
    
    plots = c(plots, list(genodepth_raw, genodepth_log))
    
}



###
# Site summary
###

if(!is.null(args$sitefile)){
    
    cat("Processing site data from",args$sitefile, "\n")
    
    sites = read.delim(args$sitefile, check.names=F)
    
    if(!is.null(args$max_datapoints) && nrow(sites) > args$max_datapoints){ # Subset if needed
        cat("\tRandomly subsetting to",args$max_datapoints,"data points\n")
        sites = sites[sample(1:nrow(sites), args$max_datapoints, replace=F), ]
    }
    
    
    # Percent Heterozygous; for some reason qplot() was not handling this properly so did full ggplot()
    hetvals = sort(sites[,'Proportion Heterozygous'])
    hetdata=data.frame(xval=1:length(hetvals), yval=hetvals)
    hets = ggplot(hetdata, mapping=aes(x=xval, y=yval)) +
        geom_point(color="darkred") +
        labs(title="Sites - Percent Heterozygous", x="Sites (sorted)", y="Fraction Het")
    
    # Percent Missing
    missing = qplot(sites[,'Proportion Missing'], bins=50, fill=I('darkgreen')) +
        labs(title="Sites - Percent Missing", x="Fraction Missing", y="# Sites")
    
    # Minor allele frequency
    maf = qplot(x=sites[,'Minor Allele Frequency'], bins=50, fill=I('darkblue')) +
        labs(title="Sites - Minor Allele Frequency", x="MAF", y="# Sites")
    

    plots = c(plots, list(hets, missing, maf))
        
}


###
# Taxa summary 
###

if(!is.null(args$taxafile)){
    
    cat("Processing taxa data from",args$sitefile, "\n")
    
    taxa = read.delim(args$taxafile, check.names=F)
    
    # Percent Heterozygous
    hetvals = sort(taxa[,'Proportion Heterozygous'])
    hets = qplot(hetvals, bins=50, fill=I('coral')) +
        labs(title = "Taxa - Percent Heterozygous", x="Fraction Het", y="# Taxa") 
        
    # Percent Missing
    missing = qplot(taxa[,'Proportion Missing'], bins=50, fill=I('aquamarine3')) +
        labs(title="Taxa - Percent Missing", x="Fraction Missing", y="# Taxa")
    
    
    
    
    plots = c(plots, list(hets, missing))
    
}


###
# Output
###

nplots = length(plots)
ncol=2
nrow = ceiling(nplots/ncol)

cat("Writing graphic to",args$outfile, "\n")
png(args$outfile, height=3*nrow, width=4*ncol, units="in", res=150)
    grid.arrange(grobs=plots, ncol=ncol, nrow=nrow, as.table=TRUE)
dev.off()
