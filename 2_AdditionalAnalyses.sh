#! /bin/bash

# This script contains additional analyses requested by reviewers during publication

TASSEL5="run_pipeline.pl -Xms10g -Xmx40g"
hplc_results="0_HDGS_HPLC_Results.csv"
genos=1g_hemp_filtered.sorted.vcf.gz

workdir="2_RevisionAdditions"
if [ ! -e $workdir ] ; then mkdir $workdir ; fi


###################
# CONDA ENVIRONMENT
###################

# Conda environment name. See accompanying 0_CreateCondaEnvironment.sh file for how to create an identical environment
conda_hemp=proj-hemp-genotyping-2020  

# Load functions required to be able to activate conda environments within scripts. KEEP THIS COMMAND UNCOMMENTED
. $(conda info --root)/etc/profile.d/conda.sh   
conda activate $conda_hemp


# ##############
# # Genetic versus Chemical diversity
# ##############

# $TASSEL5 -vcf $genos -distanceMatrix -export $workdir/2a_hemp_distances.txt
Rscript $scriptdir/2a_CompareGeneticVersusChemicalDiversity.r -c $hplc_results -d $workdir/2a_hemp_distances.txt -k $scriptdir/0_HDGS_keyfile.txt -o 2b_compare_genetics_chemistry


####################
# Add additional hemp genotypes to the dendrogram
####################

# Potential sources
# PRJNA285813 - From Sawler et al 2015(doi: 10.1371/journal.pone.0133292). GBS, but with ApeKI, so probably not much overlap
# PRJNA575581 & SUB6635057 - McKernan et al 2020 (https://doi.org/10.1101/2020.01.03.894428; bioRxiv preprint) - WGS for 3 lines and transcriptome for ~40 others
# PRJNA713792 - NCBI, sampes from Netherlands (for sale, likely); unable to find a publication
# PRJNA683613 - WGS of 23 varieties; unable to find publication? Also looks like may not be released yet
# PRJNA562042 & PRJNA561069 - Draft genome of wild Cannabis from China
# PRJNA560384 - Reference genome
# PRJNA551485 - Transcriptomes from 2 cultivars
#*PRJNA470994 - WGS(?) of "Cultivation classic" cultivars by Phylos Bioscience
#*PRJNA419020 - Diversity analysis of Iranian cannabis (https://www.nature.com/articles/s41598-017-15816-5)
#*PRJNA347566 - More Phylos biosciences
# PRJNA318063 - Reference genome again
# PRJNA317659 - GBS data of 180+ strains; unknown methodology
# PRJNA310948 - WGS of 57(?) hemp strains (Lynch et al 2017; https://doi.org/10.1080/07352689.2016.1265363)


# TODO: Find a way to doanload data, align to genome, slice out needed SNPs, and see what aligns. Want SRA data tables and ability to just pull from NCBI and align in one go.
#  Look up Kraken classification pipeline; probably can base it off that.
