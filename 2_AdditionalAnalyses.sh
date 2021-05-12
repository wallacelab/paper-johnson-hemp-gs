#! /bin/bash

# This script contains additional analyses requested by reviewers during publication
# Additional WGS of hemp taken from PRJNA310948, which includes WGS of 55 hemp varieties (Lynch et al 2017; https://doi.org/10.1080/07352689.2016.1265363)

TASSEL5="run_pipeline.pl -Xms10g -Xmx40g"

# Existing data files/directories
scriptdir=0_Scripts
hplc_results="0_HDGS_HPLC_Results.csv"
genos=1g_hemp_filtered.sorted.vcf.gz
sra_run_info=$scriptdir/PRJNA310948_run_info.csv
genome="$HOME/Projects/0_RawData/HempGenome/GCF_900626175.1_cs10_genomic.fna"
ncores=7

# Working directory for these analyses
workdir="2_RevisionAdditions"
fastqdir=$workdir/2c_fastq
if [ ! -e $workdir ] ; then mkdir $workdir ; fi
if [ ! -e $fastqdir ]; then mkdir $fastqdir; fi


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
# Rscript $scriptdir/2a_CompareGeneticVersusChemicalDiversity.r --chemistry $hplc_results --distances $workdir/2a_hemp_distances.txt --keyfile $scriptdir/0_HDGS_keyfile.txt \
#    --remove-outliers Baox -o $workdir/2a_compare_genetics_chemistry 

   
####################
# Neighbor Net instead of neighbor joinging
####################
   
# Rscript $scritpdir/2b_NieghborNetOfPlants.r --distances $workdir/2a_hemp_distances.txt --keyfile $scriptdir/0_HDGS_keyfile.txt -o $workdir/2b

####################
# Add additional hemp genotypes to the dendrogram
####################

# # For computation's sake, only get sites that made it through our genotyping. Indicate 1 bp to either side to be sure
# zcat $genos | grep -v "^#" | cut -f1-2 | sed -r "s|(.+)\t(.+)|\1\t\2|" > $workdir/2c_sites.txt
# Rscript -e "x=read.delim('$workdir/2c_sites.txt', header=F); x[,3] = x[,2] +1; x[,2] = x[,2]-1" \
#     -e "write.table(x, file='$workdir/2c_targets.bed', sep=' ', quote=F, row.names=F, col.names=F)"


# # Identify the runs to keep
# Rscript -e "runs=read.csv('$sra_run_info', stringsAsFactors=F); runs=subset(runs, runs\$LibraryStrategy=='WGS')" \
#     -e "write.table( runs[,c('Run','LibraryName')], file='$workdir/2c_sra_runs.txt', quote=T, row.names=F, col.names=F)" 

# # Index genome
# bwa index $genome

# # For each run, download, align to genome, and keep only selected sites so as to not have huge overhead
# while read run name; do
#     run=${run//\"/} # Strip quotes from accession
#     echo "Processing $run"
# 
#     # Download if not already retrieved
#     if [ ! -e $fastqdir/${run}_1.fastq ] && [ ! -e $fastqdir/${run}_1.fastq.gz ] ; then
#         echo -e "\tDownloading FASTQ files..."
#         fastq-dump --clip --split-files $run --outdir $fastqdir
#     fi
#     
#     # Align
#     bwa mem -t $ncores $genome $fastqdir/${run}_1.fastq $fastqdir/${run}_2.fastq | \
#         samtools view -F 0x04 -b -h -S -L $workdir/2c_sites.bed - > $fastqdir/$run.targets.bam
# 
#     # Gzip
#     gzip $fastqdir/${run}_1.fastq $fastqdir/${run}_2.fastq
# 
# #     break
# 
# done < $workdir/2c_sra_runs.txt


# # Sort & index bam files
# for bam in $fastqdir/*.targets.bam; do
#     basename=${bam/.targets.bam/}
#     samtools sort $bam $basename.sorted
#     samtools index $basename.sorted.bam
# done

# # Call SNPs using same options as our data
# max_depth=10000
# min_base_quality=20 
# bcftools mpileup -d $max_depth -f "$genome" -Q $min_base_quality --output-type u --skip-indels --threads $ncores --annotate FORMAT/AD,FORMAT/DP $fastqdir/*.sorted.bam  | \
#     bcftools call --multiallelic-caller --threads $ncores --skip-variants indels --output $workdir/2d_ncbi.snp_calls.bcf --output-type b --variants-only 

# # Convert to VCF
# bcftools view $workdir/2d_ncbi.snp_calls.bcf --output-type z > $workdir/2d_ncbi.snp_calls.vcf.gz

# Find common sites between the datasets
# zcat $genos | grep -v "^##" | cut -f1-5 > $workdir/2e_hdgs_sites.txt
# zcat $workdir/2d_ncbi.snp_calls.vcf.gz | grep -v "^##" | cut -f1-5 > $workdir/2e_ncbi_sites.txt
# Rscript $scriptdir/2e_CompareGenoSites.r -a $workdir/2e_hdgs_sites.txt -b  $workdir/2e_ncbi_sites.txt -o  $workdir/2e_shared_sites.txt

# Subset each dataset to shared sites
# bcftools index $workdir/2d_ncbi.snp_calls.vcf.gz
# bcftools filter --regions-file $workdir/2e_shared_sites.txt $workdir/2d_ncbi.snp_calls.vcf.gz | bgzip > $workdir/2f_ncbi.shared_snps.vcf.gz
# zcat $genos | bgzip > $workdir/2f_hdgs.all_snps.vcf.gz
# bcftools index $workdir/2f_hdgs.all_snps.vcf.gz
# bcftools filter --regions-file $workdir/2e_shared_sites.txt $workdir/2f_hdgs.all_snps.vcf.gz | bgzip > $workdir/2f_hdgs.shared_snps.vcf.gz

# # Merge; using TASSEL because is more robust that BCFtoools, but does lose a bunch of annotation data
# $TASSEL5 -fork1 -vcf $workdir/2f_hdgs.shared_snps.vcf.gz  -fork2 -vcf $workdir/2f_ncbi.shared_snps.vcf.gz -combine3 -input1 -input2 -mergeGenotypeTables \
#     -export $workdir/2g_combined_genos.vcf.gz -exportType VCF

# Stopped here because a quick manual inspection showed that the two datasets were clearly and distinctly separating from each other, likely due to the much lower depth of the NCBI data causing artifacts
