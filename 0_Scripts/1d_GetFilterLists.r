#! /usr/bin/Rscript

# Apply filters to determine which sites/taxa to keep

# Libraries
library(argparse)

# Command-line arguments
options(stringsAsFactors=F)
parser=ArgumentParser()
parser$add_argument("-s", "--sitefile", help="TASSEL genotype summary file of sites")
parser$add_argument("-t", "--taxafile", help="TASSEL genotype summary file of taxa")
parser$add_argument("-d", "--depthfile", help="VCFTools depth output file")
parser$add_argument("--site-min-depth", type="integer", default=0, help="Minimum total depth to keep a site")
parser$add_argument("--site-max-depth", type="integer", default=9999999, help="Maximum total depth to keep a site")
parser$add_argument("--site-min-het", type="double", default=0, help="Minimum fraction het to keep a site")
parser$add_argument("--site-max-het", type="double", default=1, help="Maximum fraction het to keep a site")
parser$add_argument("--site-min-maf", type="double", default=0, help="Minimum minor allele frequency to keep a site")
parser$add_argument("--site-max-maf", type="double", default=1, help="Maximum minor allele frequency to keep a site")
parser$add_argument("--site-max-missing", type="double", default=1, help="Maximum fraction missing data to keep a site")
parser$add_argument("--taxa-max-missing", type="double", default=1, help="Maximum fraction missing data to keep a taxon/sample")
parser$add_argument("--outtaxa", help="Output file of taxa to keep")
parser$add_argument("--outsites", help="Output file of sites to keep (coordinates)")
parser$add_argument("--outsitenames", help="Output file of sites to keep (names)")
args=parser$parse_args()
# setwd('/home/jgwall/Projects/PearlMillet/BoubacarPearlMilletPops/GenomeLinkageMaps/2_Analysis/')
# args=parser$parse_args(c("-s","2a_sitesummary.txt", "-t", "2a_taxasummary.txt", "-d", "2a_depth.txt",
#                          "--outtaxa", "99_tmp.taxa.txt", "--outsite", "99_tmp.site.txt",
#                          "--site-max-depth", "1000", "--site-max-het", "0.25", "--site-min-maf", "0.2",
#                          "--taxa-max-missing", "0.9", "--site-max-missing", "0.6"))


# Load data
depths = read.delim(args$depthfile, check.names=F)
sites = read.delim(args$sitefile, check.names=F)
taxa = read.delim(args$taxafile, check.names=F)

#### 
# Filter sites
####

cat("Compiling site filters on", nrow(sites), "sites\n")

# Depth
# if(!identical(depths$ID, sites$`Site Name`)){ # Old check for site names; bcftools doesn't generate them like TASSEL, though, so doesn't work.
#     stop("Depth info and site info names do not match\n")
# }
good_depth = depths$SUM_DEPTH >= args$site_min_depth & depths$SUM_DEPTH <= args$site_max_depth
cat("\t",sum(good_depth, na.rm=T),"sites pass depth filter of [", args$site_min_depth, ",", args$site_max_depth, "]\n")

# Het
good_het = sites$`Proportion Heterozygous` >= args$site_min_het & sites$`Proportion Heterozygous` <= args$site_max_het
cat("\t",sum(good_het, na.rm=T),"sites pass het filter of [", args$site_min_het, ",", args$site_max_het, "]\n")

# MAF
good_maf = sites$`Minor Allele Frequency` >= args$site_min_maf & sites$`Minor Allele Frequency` <= args$site_max_maf
cat("\t",sum(good_maf, na.rm=T),"sites pass MAF filter of [", args$site_min_maf, ",", args$site_max_maf, "]\n")

# Missing
good_missing = sites$`Proportion Missing` <= args$site_max_missing
cat("\t",sum(good_missing, na.rm=T),"sites pass missingness filter of [ 0 , ", args$site_max_missing, "]\n")

# Combine
good_sites = good_depth & good_het & good_maf & good_missing
cat("\t",sum(good_sites),"sites are left after combining all filters\n")

# Output
output = sites[good_sites, c('Chromosome', 'Physical Position')]
# # Confirm uniqueness of each chrom-pos
chrpos = paste(sites$Chromosome, sites$`Physical Position`)
if(length(table(table(chrpos))) != 1){
    cat("\t\tWARNING!! Chromosome-position combinations do not look unique!")
}
cat("\tOutputting good sites as chrom-position to", args$outsites, "\n")
write.table(output, file=args$outsites, sep='\t', quote=F, row.names=F, col.names=F)

# Optional output of names for sites, not just coordinates
if(!is.null(args$outsitenames)){
    cat("\tOutputting good sites' names to", args$outsitenames, "\n")
    sitenames = sites[good_sites, 'Site Name']
    write(sitenames, file=args$outsitenames)
}


#### 
# Filter taxa
####

cat("Compiling taxa filters on", nrow(taxa), "taxa\n")

# Missing
good_taxmissing = taxa$`Taxa Name`[taxa$`Proportion Missing` <= args$taxa_max_missing]
cat("\t",length(good_taxmissing),"taxa pass missingness filter of [ 0 , ", args$taxa_max_missing, "]\n")

# Combine [kept here in case expand script later]
good_taxa = good_taxmissing
cat("\t",length(good_taxa),"taxa are left after combining all filters\n")

cat("\tOutputting good taxa names to", args$outtaxa, "\n")
write(good_taxa, file=args$outtaxa, sep='\n')
