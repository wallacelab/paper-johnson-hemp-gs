#! /usr/bin/Rscript

# Compare hemp genetics and chemistry to see if more variable in one is more vairable in other

# Libraries
library(argparse)

# Command-line arguments
options(stringsAsFactors=F)
parser=ArgumentParser()
parser$add_argument("-a", "--a", help="VCF site data (first 5 columns) of file A")
parser$add_argument("-b", "--b", help="VCF site data (first 5 columns) of file B")
parser$add_argument("-o", "--outfile", help="Output file of common sites (names only)")
args=parser$parse_args()
# setwd('/home/jgwall/Projects/Hemp/HempAnalysis_JasonRedo/2_RevisionAdditions/')
# args=parser$parse_args(c("-a", "2e_hdgs_sites.txt", "-b", "2e_ncbi_sites.txt", "-o", "99_tmp.txt"))

# Load data
cat("Loading VCF site info:\n\tA: ",args$a,"\n\tB:",args$b, "\n")
a = read.delim(args$a)
b = read.delim(args$b)
names(a)[1] = names(b)[1] = "CHROM" # correct reading of '#'

# Compare chrom + position
a$site = paste(a$CHROM, a$POS, sep=":")
b$site = paste(b$CHROM, b$POS, sep=":")
common=intersect(a$site, b$site)
cat("\tSet A has",nrow(a),"sites\n")
cat("\tSet B has",nrow(b),"sites\n")
cat("\tThere are",length(common),"sites in common\n")

# Subset to just shared sites and sort
shared.a = subset(a, a$site %in% common)
shared.b = subset(b, b$site %in% common)
shared.a = shared.a[order(shared.a$CHROM, shared.a$POS),]
shared.b = shared.b[order(shared.b$CHROM, shared.b$POS),]

# Compare alleles
ref_match   = sum(shared.a$REF == shared.b$REF)
ref_nomatch = sum(shared.a$REF != shared.b$REF)
alt_match   = sum(shared.a$ALT == shared.b$ALT)
alt_nomatch = sum(shared.a$ALT != shared.b$ALT)
cat("\t",ref_match,"REF alleles match and",ref_nomatch,"do not match\n")
cat("\t",alt_match,"ALT alleles match and",alt_nomatch,"do not match\n")

# Keep only the sites where the alternate alleles match
cat("\tKeeping only",alt_match,"sites where alternate alleles match\n")
mismatch = shared.a$ALT != shared.b$ALT
tokeep = shared.a[!mismatch,c("CHROM","POS")]

# Write out shared sites for future filtering
write.table(tokeep, file=args$outfile, sep="\t", quote=F, row.names=F, col.names=F)
