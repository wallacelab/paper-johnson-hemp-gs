#! /bin/bash

# Rerunning the hemp genotyping pipeline because first run filtered out PCR duplicates, which is all GBS is

# At that depth, het calls from sequencing errors will be common. Need to find the best way to call those. Can GATK do this well?
#    I wonder if I'll need to code my own het caller. Graph % minor allele for all het calls

# TODO FIRST: Check het calls to see what fraction of minor allele needed to make a call (when depth is high?)  -> Seem okay. Average around 50%, min call of each of 3. May be more aggressive.

# TODO: Correct het calls based on simple probabilities of allele depths, not HWE ?  -> See how they come out, then correct if need be

TASSEL5="run_pipeline.pl -Xms10g -Xmx40g"
bam_dir="/media/jgwall/Seagate Expansion Drive/Hemp_diversity_GS/Hemp_variant_calling/Bams/"
genome="$HOME/Projects/0_RawData/HempGenome/GCF_900626175.1_cs10_genomic.fna"
restriction_site="CAYNNNNRTG" # MslI restriction site in normal nucleotides. Mostly here for reference
hplc_results="0_HDGS_HPLC_Results.csv"
ncores=7 # number of processor cores to use

# Subdirectories
scriptdir="0_Scripts"
sortdir="1_bams_sorted"
filtdir="1a_filtering"
if [ ! -e $sortdir ] ; then mkdir $sortdir; fi
if [ ! -e $filtdir ] ; then mkdir $filtdir; fi

###################
# CONDA ENVIRONMENT
###################

# Conda environment name. See accompanying 0_CreateCondaEnvironment.sh file for how to create an identical environment
conda_hemp=proj-hemp-genotyping-2020  

# Load functions required to be able to activate conda environments within scripts. KEEP THIS COMMAND UNCOMMENTED
. $(conda info --root)/etc/profile.d/conda.sh   
conda activate $conda_hemp

##############
# Rebuild genotypes from BAM files
##############

# for sam in "$bam_dir"/*.aligned.sam; do
#     name=`basename "$sam"`
#     name=${name/.sam/}
#     echo $name
# 
#       
#     samtools view -h -S -b "$sam" > $sortdir/$name.bam
#     samtools sort -@ $ncores $sortdir/$name.bam $sortdir/$name.sorted
#     samtools index $sortdir/$name.sorted.bam
#     rm $sortdir/$name.bam
#     
# done

# # Call SNPs on all samples (skip indels with --skip-variants b/c these are harder to score and can cause issues). This takes ~36 hours on an 8-core machine
# max_depth=10000
# min_base_quality=20 
# bcftools mpileup -d $max_depth -f "$genome" -Q $min_base_quality --output-type u --skip-indels --threads $ncores --annotate FORMAT/AD,FORMAT/DP $sortdir/*.sorted.bam  | \
#     bcftools call --multiallelic-caller --threads $ncores --skip-variants indels --output 1b_hemp.snp_calls.bcf --output-type b --variants-only 

##############
# Basic Filtering - Depth
##############

# bcftools view 1b_hemp.snp_calls.bcf --output-type z > 1b_hemp.snp_calls.vcf.gz
# 
# # Get reports
# $TASSEL5 -SortGenotypeFilePlugin -inputFile 1b_hemp.snp_calls.vcf.gz -outputFile 1b_hemp.snp_calls.sorted.vcf.gz
# $TASSEL5 -vcf 1b_hemp.snp_calls.sorted.vcf.gz -genotypeSummary site -export $filtdir/1c_sitesummary.txt
# $TASSEL5 -vcf 1b_hemp.snp_calls.sorted.vcf.gz -genotypeSummary taxa -export $filtdir/1c_taxasummary.txt
# vcftools --gzvcf 1b_hemp.snp_calls.sorted.vcf.gz  --site-depth --stdout > $filtdir/1c_site_depth.txt

# # Complicated bash one-liner to count size of predicted restriction fragments in the hemp genome
# cat "$genome" | sed -r -e "s|^>.+|=====|" | tr -d '\n' | tr [:lower:] [:upper:] | \
#     sed -r -e "s|(CA[CT]..)(..[AG]TG)|\1\n\2|g" -e "s|=====||g" | \
#     awk '{print length}' > $filtdir/1c_restriction_frag_lengths.txt

# # Plot diagnostic plots from the above summary data
# Rscript $scriptdir/1c_PlotGenoSummary.r --sitefile $filtdir/1c_sitesummary.txt --taxafile $filtdir/1c_taxasummary.txt --depthfile $filtdir/1c_site_depth.txt \
#     --restriction-frags $filtdir/1c_restriction_frag_lengths.txt --max-datapoints 10000 -o $filtdir/1c_summary_plots.png


# Hmm...based on number of restriction fragments and read depth, I'd expect a true mean depth of single-copy DNA around 5, but actual cluster seems around 100. 
#   (Ignoring the cluster of aligment errors that only have 1 read) Let's filter for the 100x cluster and see how it works


# # Initial filter for mean depth and % presence
# site_min_mean_depth=15
# site_max_mean_depth=125
# site_max_missing=0.2  # Equivalent to 80% present
# Rscript $scriptdir/1d_GetFilterLists.r --sitefile $filtdir/1c_sitesummary.txt --taxafile $filtdir/1c_taxasummary.txt --depthfile $filtdir/1c_site_depth.txt \
#     --site-min-mean-depth $site_min_mean_depth --site-max-mean-depth $site_max_mean_depth --site-max-missing $site_max_missing \
#     --outtaxa $filtdir/1d_taxa_to_keep.txt --outsites $filtdir/1d_sites_to_keep.txt

# # Perform actual filtering
# bcftools view --targets-file $filtdir/1d_sites_to_keep.txt --samples-file $filtdir/1d_taxa_to_keep.txt --output-type z --output-file $filtdir/1e_hemp_filter1.vcf.gz 1b_hemp.snp_calls.bcf 


################
# Check (and correct) het calls
################

# # Heterozygous calls are a minor part of the genotypes, but I want to do a quick check to make sure they're reasonable.
# # Note: This code as written removes tertiary allele states. They affect some genotypes, but not many
# echo -e "call\tref\talt" > $filtdir/1f_check_het_calls.txt
# zcat $filtdir/1e_hemp_filter1.vcf.gz | grep -v "^#" | cut -f10- | tr '\t' '\n' | tr ":" "\t", | cut -f1-2 | grep -e "0/1" -e "1/0" | tr ',' '\t' | cut -f1-3 >> $filtdir/1f_check_het_calls.txt
# Rscript $scriptdir/1f_CheckHetCalls.r -i $filtdir/1f_check_het_calls.txt -o $filtdir/1f_check_het_calls.png --log-scale

# Okay, definitely looks like the genotype caller is being too generous with het calls, probably b/c it assumes Hardy-Weinberg Equilibrium. Recall hets (and only hets) based on binomial probabilities
# python3 $scriptdir/1g_CorrectHetsInVcf.py -i $filtdir/1e_hemp_filter1.vcf.gz --p-cutoff 0.01 -o $filtdir/1g_hemp_filter1.fix_hets.vcf.gz #--debug


##############
# Final Filtering - Allele frequency and heterozygosity
##############


# # Run post-filter diagnostic plots (really should make this step its own function or subscript or something)
# $TASSEL5 -SortGenotypeFilePlugin -inputFile $filtdir/1g_hemp_filter1.fix_hets.vcf.gz -outputFile $filtdir/1g_hemp_filter1.fix_hets.sorted.vcf.gz
# $TASSEL5 -vcf $filtdir/1g_hemp_filter1.fix_hets.sorted.vcf.gz -genotypeSummary site -export $filtdir/1h_sitesummary.filter1.txt
# $TASSEL5 -vcf $filtdir/1g_hemp_filter1.fix_hets.sorted.vcf.gz -genotypeSummary taxa -export $filtdir/1h_taxasummary.filter1.txt
# vcftools --gzvcf $filtdir/1g_hemp_filter1.fix_hets.sorted.vcf.gz  --site-depth --stdout > $filtdir/1h_site_depth.filter1.txt
# 
# # Plot 
# Rscript $scriptdir/1c_PlotGenoSummary.r --sitefile $filtdir/1h_sitesummary.filter1.txt --taxafile $filtdir/1h_taxasummary.filter1.txt --depthfile $filtdir/1h_site_depth.filter1.txt \
#    --restriction-frags $filtdir/1c_restriction_frag_lengths.txt --max-datapoints 10000 -o $filtdir/1h_hemp_filter1.summary_plots.png


    
# # Second round of filtering
# site_max_het=0.1
# site_min_maf=0.025
# Rscript $scriptdir/1d_GetFilterLists.r --sitefile $filtdir/1h_sitesummary.filter1.txt --taxafile $filtdir/1h_taxasummary.filter1.txt --depthfile $filtdir/1h_site_depth.filter1.txt \
#     --site-max-het $site_max_het --site-min-maf $site_min_maf \
#     --outtaxa $filtdir/1i_taxa_to_keep.filt2.txt --outsites $filtdir/1i_sites_to_keep.filt2.txt

# # Perform actual filtering
# bcftools view --targets-file $filtdir/1i_sites_to_keep.filt2.txt --samples-file $filtdir/1i_taxa_to_keep.filt2.txt --output-type z --output-file $filtdir/1j_hemp_filter2.vcf.gz $filtdir/1g_hemp_filter1.fix_hets.vcf.gz


# # # Final post-filter diagnostic plots (mostly as a sanity check)
# $TASSEL5 -SortGenotypeFilePlugin -inputFile $filtdir/1j_hemp_filter2.vcf.gz -outputFile $filtdir/1j_hemp_filter2.sorted.vcf.gz
# $TASSEL5 -vcf $filtdir/1j_hemp_filter2.sorted.vcf.gz -genotypeSummary site -export $filtdir/1j_sitesummary.filter2.txt
# $TASSEL5 -vcf $filtdir/1j_hemp_filter2.sorted.vcf.gz -genotypeSummary taxa -export $filtdir/1j_taxasummary.filter2.txt
# vcftools --gzvcf $filtdir/1j_hemp_filter2.sorted.vcf.gz  --site-depth --stdout > $filtdir/1j_site_depth.filter2.txt
# 
# # Plot 
# Rscript $scriptdir/1c_PlotGenoSummary.r --sitefile $filtdir/1j_sitesummary.filter2.txt --taxafile $filtdir/1j_taxasummary.filter2.txt --depthfile $filtdir/1j_site_depth.filter2.txt \
#     --restriction-frags $filtdir/1c_restriction_frag_lengths.txt --max-datapoints 10000 -o $filtdir/1j_hemp_filter2.summary_plots.png



# Copy filtered genotypes back to main directory
cp $filtdir/1j_hemp_filter2.sorted.vcf.gz ./1g_hemp_filtered.sorted.vcf.gz



##############
# Phylogenetic tree
##############

# # Tree
# $TASSEL5 -vcf 1g_hemp_filtered.sorted.vcf.gz -tree Neighbor -treeSaveDistance false -export 1g_hemp_filtered.sorted.tre -exportType Text
# mv 1g_hemp_filtered.sorted.tre.txt 1g_hemp_filtered.sorted.tre


##############
# Genomic prediction
##############

# # Rename samples to match phenotype info (Tried with 'bcftools reheader', but it did the compression wrong somehow)
source=1g_hemp_filtered.sorted.vcf.gz
newfile=1h_hemp_filtered.sorted.renamed.vcf
zcat $source | grep "^##" > $newfile   # Metadata lines
zcat $source | grep "^#CHROM" | sed -r "s|1_bams_sorted/P.-HDGS-([0-9]+).aligned.sorted.bam|HDGS_\1|g" >> $newfile  # Header line
zcat $source | grep -v "^#" >> $newfile   # Genotypes
gzip $newfile

# Put phenotypes in TASSEL format; also do sqrt and log transformations, just in case
Rscript $scriptdir/1i_FormatPhenosForTassel.r -i $hplc_results -o 1i_hdgs_phenos --log --sqrt

# Run genomic prediction with TASSEL
$TASSEL5 -fork1 -vcf 1h_hemp_filtered.sorted.renamed.vcf.gz -ck -fork2 -t 1i_hdgs_phenos.tassel.txt -combine3 -input1 -input2 \
    -GenomicSelectionPlugin -doCV true -kFolds 10 -nIter 1000 -endPlugin -export 1j_gs_results.txt
    
# Graph results as violin plots
Rscript $scriptdir/1k_PlotPredictionAccuracy.r -i 1j_gs_results.txt -o 1k_gs_results

# Not great, unfortunately



# Below for backup purposes
# # # # # # ##############
# # # # # # # Final Filtering - Allele frequency and heterozygosity
# # # # # # ##############
# # # # # # 
# # # # # # 
# # # # # # 
# # # # # # 
# # # # # # # # Run post-filter diagnostic plots (really should make this step its own function or subscript or something)
# # # # # # # $TASSEL5 -SortGenotypeFilePlugin -inputFile $filtdir/1e_hemp_filter1.vcf.gz -outputFile $filtdir/1e_hemp_filter1.sorted.vcf.gz
# # # # # # # $TASSEL5 -vcf $filtdir/1e_hemp_filter1.sorted.vcf.gz -genotypeSummary site -export $filtdir/1e_sitesummary.filter1.txt
# # # # # # # $TASSEL5 -vcf $filtdir/1e_hemp_filter1.sorted.vcf.gz -genotypeSummary taxa -export $filtdir/1e_taxasummary.filter1.txt
# # # # # # # vcftools --gzvcf $filtdir/1e_hemp_filter1.sorted.vcf.gz  --site-depth --stdout > $filtdir/1e_site_depth.filter1.txt
# # # # # # 
# # # # # # # # Plot 
# # # # # # # Rscript $scriptdir/1c_PlotGenoSummary.r --sitefile $filtdir/1e_sitesummary.filter1.txt --taxafile $filtdir/1e_taxasummary.filter1.txt --depthfile $filtdir/1e_site_depth.filter1.txt \
# # # # # # #    --restriction-frags $filtdir/1c_restriction_frag_lengths.txt --max-datapoints 10000 -o $filtdir/1e_hemp_filter1.summary_plots.png
# # # # # # 
# # # # # # 
# # # # # #     
# # # # # # # # # Second round of filtering
# # # # # # # site_max_het=0.1
# # # # # # # site_min_maf=0.025
# # # # # # # Rscript $scriptdir/1d_GetFilterLists.r --sitefile $filtdir/1e_sitesummary.filter1.txt --taxafile $filtdir/1e_taxasummary.filter1.txt --depthfile $filtdir/1e_site_depth.filter1.txt \
# # # # # # #     --site-max-het $site_max_het --site-min-maf $site_min_maf \
# # # # # # #     --outtaxa $filtdir/1f_taxa_to_keep.filt2.txt --outsites $filtdir/1f_sites_to_keep.filt2.txt
# # # # # # # 
# # # # # # # # Perform actual filtering
# # # # # # # bcftools view --targets-file $filtdir/1f_sites_to_keep.filt2.txt --samples-file $filtdir/1f_taxa_to_keep.filt2.txt --output-type z --output-file $filtdir/1g_hemp_filter2.vcf.gz $filtdir/1e_hemp_filter1.vcf.gz
# # # # # # 
# # # # # # 
# # # # # # # # # Final post-filter diagnostic plots (mostly as a sanity check)
# # # # # # # $TASSEL5 -SortGenotypeFilePlugin -inputFile $filtdir/1g_hemp_filter2.vcf.gz -outputFile $filtdir/1g_hemp_filter2.sorted.vcf.gz
# # # # # # # $TASSEL5 -vcf $filtdir/1g_hemp_filter2.sorted.vcf.gz -genotypeSummary site -export $filtdir/1g_sitesummary.filter2.txt
# # # # # # # $TASSEL5 -vcf $filtdir/1g_hemp_filter2.sorted.vcf.gz -genotypeSummary taxa -export $filtdir/1g_taxasummary.filter2.txt
# # # # # # # vcftools --gzvcf $filtdir/1g_hemp_filter2.sorted.vcf.gz  --site-depth --stdout > $filtdir/1g_site_depth.filter2.txt
# # # # # # 
# # # # # # # # Plot 
# # # # # # # Rscript $scriptdir/1c_PlotGenoSummary.r --sitefile $filtdir/1g_sitesummary.filter2.txt --taxafile $filtdir/1g_taxasummary.filter2.txt --depthfile $filtdir/1g_site_depth.filter2.txt \
# # # # # # #     --restriction-frags $filtdir/1c_restriction_frag_lengths.txt --max-datapoints 10000 -o $filtdir/1g_hemp_filter2.summary_plots.png
# # # # # # 
# # # # # # 
# # # # # # 
# # # # # # # Copy filtered genotypes back to main directory
# # # # # # # cp $filtdir/1g_hemp_filter2.sorted.vcf.gz ./1g_hemp_filtered.sorted.vcf.gz
