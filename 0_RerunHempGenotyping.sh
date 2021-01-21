#! /bin/bash

# Rerunning the hemp genotyping pipeline because first run filtered out PCR duplicates, which is all GBS is

# TODO: Looking at mean depth vs coverage, I _Think_ an average coverage of 100x is about right. However, need to confirm.
# At that depth, het calls from sequencing errors will be common. Need to find the best way to call those. Can GATK do this well?
#    I wonder if I'll need to code my own het caller. Graph % minor allele for all het calls

# Command to pull out all het genotype calls: zcat 1b_hemp.snp_calls.sorted.vcf.gz | grep -v "^#" | cut -f10-12 | tr '\t' '\n' | tr ":" "\t", | cut -f1-2 | grep -e "0/1" -e "1/0" > tmp_hets.txt
#   The above is useful to ouble check calls. To get just ref/alt, add  "cut -f2 | tr ',' '\t'"
# TODO FIRST: Check het calls to see what fraction of minor allele needed to make a call (when depth is high?)  -> Seem okay. Average around 50%, min call of each of 3. May be more aggressive.

TASSEL5="run_pipeline.pl -Xms10g -Xmx40g"
bam_dir="/media/jgwall/Seagate Expansion Drive/Hemp_diversity_GS/Hemp_variant_calling/Bams/"
genome="/media/jgwall/Seagate Expansion Drive/Hemp_diversity_GS/Hemp_variant_calling/GCF_900626175.1_cs10_genomic.fna"
keyfile=0_HDGS_SplitsTree_KeyFile.txt
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
# Basic Filtering
##############

# bcftools view 1b_hemp.snp_calls.bcf --output-type z > 1b_hemp.snp_calls.vcf.gz
# 
# # Get reports
# $TASSEL5 -SortGenotypeFilePlugin -inputFile 1b_hemp.snp_calls.vcf.gz -outputFile 1b_hemp.snp_calls.sorted.vcf.gz
# $TASSEL5 -vcf 1b_hemp.snp_calls.sorted.vcf.gz -genotypeSummary site -export $filtdir/1c_sitesummary.txt
# $TASSEL5 -vcf 1b_hemp.snp_calls.sorted.vcf.gz -genotypeSummary taxa -export $filtdir/1c_taxasummary.txt
# vcftools --gzvcf 1b_hemp.snp_calls.sorted.vcf.gz  --site-depth --stdout > $filtdir/1c_site_depth.txt

# # Plot diagnostic plots from the above summary data - TODO: Add depth binplot to identify sweet spot
Rscript $scriptdir/1c_PlotGenoSummary.r --sitefile $filtdir/1c_sitesummary.txt --taxafile $filtdir/1c_taxasummary.txt --depthfile $filtdir/1c_site_depth.txt \
    --max-datapoints 10000 -o $filtdir/1c_summary_plots.png

# # Basic site and individual filtering
# site_max_depth=500
# site_max_het=0.1
# site_min_maf=0.025
# site_max_missing=0.4
# Rscript $scriptdir/1d_GetFilterLists.r --sitefile $filtdir/1c_sitesummary.txt --taxafile $filtdir/1c_taxasummary.txt --depthfile $filtdir/1c_site_depth.txt \
#     --site-max-depth $site_max_depth --site-max-het $site_max_het --site-min-maf $site_min_maf \
#     --site-max-missing $site_max_missing --outtaxa $filtdir/1d_taxa_to_keep.txt --outsites $filtdir/1d_sites_to_keep.txt

# # Perform actual filtering
# bcftools view --targets-file $filtdir/1d_sites_to_keep.txt --samples-file $filtdir/1d_taxa_to_keep.txt --output-type z --output-file 1e_hemp_filtered.vcf.gz 1b_hemp.snp_calls.bcf 

# # Run post-filter diagnostic plots
# $TASSEL5 -SortGenotypeFilePlugin -inputFile 1e_hemp_filtered.vcf.gz -outputFile 1e_hemp_filtered.sorted.vcf.gz
# $TASSEL5 -vcf 1e_hemp_filtered.sorted.vcf.gz -genotypeSummary site -export $filtdir/1e_sitesummary.postfilter.txt
# $TASSEL5 -vcf 1e_hemp_filtered.sorted.vcf.gz -genotypeSummary taxa -export $filtdir/1e_taxasummary.postfilter.txt
# vcftools --gzvcf 1e_hemp_filtered.sorted.vcf.gz  --site-depth --stdout > $filtdir/1e_site_depth.postfilter.txt

# # Plot 
# Rscript $scriptdir/1c_PlotGenoSummary.r --sitefile $filtdir/1e_sitesummary.postfilter.txt --taxafile $filtdir/1e_taxasummary.postfilter.txt --depthfile $filtdir/1e_site_depth.postfilter.txt \
#    --max-datapoints 10000 -o 1e_hemp_filtered.summary_plots.png










##############
# Get distances and trees (= quick check to see how consistent these are with Matt's original calls)
##############

# orig=../HempAnalysis_MattJohnson/HDGS_Analysis/Population_Structure/HDGS_less_het_no_Indel.hmp.txt
# new=1e_hemp_filtered.sorted.vcf.gz

# Tree
# $TASSEL5 -vcf $new -tree Neighbor -treeSaveDistance false -export 1f_hemp_filtered.sorted.tre -exportType Text
# mv 1f_hemp_filtered.sorted.tre.txt 1f_hemp_filtered.sorted.tre

# # Distances
# $TASSEL5 -vcf $new -distanceMatrix -export 1f_hemp_filtered.sorted.distances.txt

# # Genotypes to hapmap format (for easier manipulation)
# $TASSEL5 -vcf $new -export 1f_hemp_filtered_.sorted.hmp.txt.gz

# # Now same for Matt's original dataset
# $TASSEL5 -h $orig -tree Neighbor -treeSaveDistance false -export 1g_hemp_orig.tre -exportType Text
# mv 1g_hemp_orig.tre.txt 1g_hemp_orig.tre
# $TASSEL5 -h $orig -distanceMatrix -export 1g_hemp_orig.distances.txt

# # Compare the two
# Rscript $scriptdir/1h_CompareHempGenos.r --oldtree 1g_hemp_orig.tre --newtree 1f_hemp_filtered.sorted.tre --olddist 1g_hemp_orig.distances.txt --newdist 1f_hemp_filtered.sorted.distances.txt \
#     --oldgenos $orig --newgenos 1f_hemp_filtered_.sorted.hmp.txt.gz  --key $keyfile --outprefix 1h_compare_old_new_genos



# Hmmm....the above is actually really bad. Only 3 sites in common, low genetic distance correlation...it's like an entirely different dataset. I should check the raw genotypes. 
# Doing raw distances in a subfolder to keep things clean


# orig=../HempAnalysis_MattJohnson/Genotyping/merged_hmps/HDGS_Merged.hmp.txt
# new=1b_hemp.snp_calls.sorted.vcf.gz

# # New Genos (raw)
# $TASSEL5 -vcf $new -tree Neighbor -treeSaveDistance false -export $rawdir/1h_newgenos_raw.tre -exportType Text
# mv $rawdir/1h_newgenos_raw.tre.txt $rawdir/1h_newgenos_raw.tre
# $TASSEL5 -vcf $new -distanceMatrix -export $rawdir/1h_newgenos_raw.distances.txt
# $TASSEL5 -vcf $new -export $rawdir/1h_newgenos_raw.hmp.txt.gz

# # Old Genos (raw)
# $TASSEL5 -h $orig -tree Neighbor -treeSaveDistance false -export $rawdir/1h_oldgenos_raw.tre -exportType Text
# mv $rawdir/1h_oldgenos_raw.tre.txt $rawdir/1h_oldgenos_raw.tre
# $TASSEL5 -h $orig -distanceMatrix -export $rawdir/1h_oldgenos_raw.distances.txt

# # Compare the two
# Rscript $scriptdir/1h_CompareHempGenos.r --oldtree $rawdir/1h_oldgenos_raw.tre --newtree $rawdir/1h_newgenos_raw.tre --olddist $rawdir/1h_oldgenos_raw.distances.txt --newdist $rawdir/1h_newgenos_raw.distances.txt \
#     --oldgenos $orig --newgenos $rawdir/1h_newgenos_raw.hmp.txt.gz  --key $keyfile --outprefix $rawdir/1h_compare_raw_genos



##############
# Redo genomic prediction to see if it's any better. (Real indicator of if these genotypes are better or not)
##############

# # Rename samples to match phenotype info (Tried with 'bcftools reheader', but it did the compression wrong somehow)
# newfile=1i_hemp_filtered.sorted.renamed.vcf
# zcat 1e_hemp_filtered.sorted.vcf.gz | grep "^##" > $newfile   # Metadata lines
# zcat 1e_hemp_filtered.sorted.vcf.gz | grep "^#CHROM" | sed -r "s|1_bams_sorted/P.-HDGS-([0-9]+).aligned.sorted.bam|HDGS_\1|g" >> $newfile  # Header line
# zcat 1e_hemp_filtered.sorted.vcf.gz | grep -v "^#" >> $newfile   # Genotypes
# gzip $newfile

# # Put phenotypes in TASSEL format; also do sqrt and log transformations, just in case
# Rscript $scriptdir/1j_FormatPhenosForTassel.r -i 0_HDGS_HPLC_Results.csv -o 1j_hdgs_phenos --log --sqrt

# # Run genomic prediction with TASSEL
# $TASSEL5 -fork1 -vcf 1i_hemp_filtered.sorted.renamed.vcf.gz -ck -fork2 -t 1j_hdgs_phenos.tassel.txt -combine3 -input1 -input2 \
#     -GenomicSelectionPlugin -doCV true -kFolds 10 -nIter 1000 -endPlugin -export 1k_gs_results.txt
    
# # Graph results as violin plots
# Rscript $scriptdir/1l_PlotPredictionAccuracy.r -i 1k_gs_results.txt -o 1l_gs_results

# Not great, unfortunately



##############
# I don't like how different the genotyping sets are. Time to look in more detail
##############

# # Check my genotypes for how consistent clones are, since that's one thing that triggered this whole investigation
# Rscript $scriptdir/1m_RenameSamples.r -i 1f_hemp_filtered_.sorted.hmp.txt.gz -k $keyfile -o 1m_hemp_filtered.sorted.renamed.hmp.txt
# # Manual inspection of the above shows mostly consistency in the BaoxSP and LifterSP clones, with the occasional het or errant genotype. Seems more consistent than Matt's original calls, though

# # Combine the raw datasets of the two and see how they compare
# orig=../HempAnalysis_MattJohnson/Genotyping/merged_hmps/HDGS_Merged.hmp.txt
# new=$rawdir/1h_newgenos_raw.hmp.txt.gz
# ## Get the sites in common between the two
# cut -f1 $orig | tail -n +2 | gzip > 1n_orig_sites.raw.txt.gz
# zcat $new | cut -f1 | tail -n +2 | gzip > 1n_new_sites.raw.txt.gz
# Rscript -e "a=scan('1n_orig_sites.raw.txt.gz', what=character()); b=scan('1n_new_sites.raw.txt.gz', what=character())" \
#     -e "both = intersect(a, b); write(both, file='1n_common_sites.raw.txt')"
# ## Filter and combine genotypes
# $TASSEL5 -fork1 -h $orig -includeSiteNamesInFile 1n_common_sites.raw.txt \
#     -fork2 -h $new -includeSiteNamesInFile 1n_common_sites.raw.txt \
#     -combine3 -input1 -input2 -mergeGenotypeTables -export 1n_orig_new.common_sites.hmp.txt.gz


# # Okay, the above is kind of disturbing when manually inspected. Not only do my new genotypes have much more coverage, in some cases the calls are directly opposite. need to check the BAM files directly.
# region1=NC_044375.1:7639468-7639511  # HDGS 100 different here (along with many others)
# region2=NC_044375.1:35153569-35153600  # HDGS 10 and 100 have different calls here (along with many others)
# p100="$bam_dir/P1-HDGS-100.aligned.bam"
# samtools sort -@ $ncores "$p100" ./1o_HDGS100.sorted
# samtools index 1o_HDGS100.sorted.bam
#samtools view -h 1o_HDGS100.sorted.bam $region1 $region2 > 1o_HDGS100.sorted.zoom.sam

# manual inspection of the above in IGVtools shows some serious problems with the old genotypes. In every case, I have an alternate allele while Matt has the reference, and the BAM files support my call
#  In one instance (NC_044375.1:35153600), Matt shows a missing call despite a read depth of 127; what's going on?
# Also, how am I getting read depths this high? I thought we were at <1x sequencing in general. Are these repetative regions?

# # Check the ones I think Matt used for his genotyping
# samtools view -h "$bam_dir/P1-HDGS-100.markedDups.fixedRG.bam" $region1 $region2 > 1o_HDGS100.markedDups.zoom.sam
# # Okay, this is confusing. There's just 1 read (rest filtered out?) but it matches my genotypes, not the ones in Matt's VCF file. how?

# # Matt used GATK; let's see if I can recreate his calls on those BAM files, and see what happens if I change the source to the non-PCR-dup-marked ones
# refGenome=/home/jgwall/Projects/0_RawData/HempGenome/GCF_900626175.1_cs10_genomic.fna

## Matt's original call, with targeted BAM, using GATK 4.1.5.0
# samtools view -bS 1o_HDGS100.markedDups.zoom.sam > 1o_HDGS100.markedDups.zoom.bam
# samtools index 1o_HDGS100.markedDups.zoom.bam
# gatk --java-options "-Xmx4g" HaplotypeCaller -R $refGenome -I 1o_HDGS100.markedDups.zoom.bam -O 1p_HDGS100_gatk_calls.markedDups.zoom.vcf -ERC GVCF -G StandardAnnotation -G AS_StandardAnnotation -G StandardHCAnnotation

## Redo with just the aligned reads, not the marked dup ones. Have to add a read group (@RG) flag for GATK to process properly
# samtools view -bS 1o_HDGS100.sorted.zoom.sam > 1o_HDGS100.sorted.zoom.bam
# picard AddOrReplaceReadGroups \
#     I=1o_HDGS100.sorted.zoom.bam \
#     O=1o_HDGS100.sorted.zoom.fix_readgroup.bam \
#     SORT_ORDER=coordinate \
#     RGID=HDGS100 \
#     RGLB=HDGS100 \
#     RGPL=illumina \
#     RGSM=HDGS100 \
#     RGPU=unk \
#     CREATE_INDEX=True
# gatk --java-options "-Xmx4g" HaplotypeCaller -R $refGenome -I 1o_HDGS100.sorted.zoom.fix_readgroup.bam -O 1p_HDGS100_gatk_calls.aligned_sorted.zoom.vcf -ERC GVCF -G StandardAnnotation -G AS_StandardAnnotation -G StandardHCAnnotation

# Okay, using the raw alignment "fixed" at least one of the regions (NC_044375.1:7639468-7639511 ) so that it shows consistent alternative alleles calls. Looks like the "mark duplicates" did indeed kill this dataset.
#   ...which means I now have to make sure _everything_ is run with the new dataset, and still double-check it to make sure things are what I think. *sigh*
# TODO: Rerun the 1o_ troubleshooting section after deleting all the files to clean it up folder's getting cluttered. Maybe put intermediate files in a subfolder
