#! /bin/bash

# BWA commands were run separately from the rest of the analysis


# specify directories and full path to ref genome for analysis
workDir=/home/mjohnson/Desktop/Hemp_diversity_GS/Hemp_variant_calling
dataDir=/home/mjohnson/Desktop/Hemp_diversity_GS/Hemp_fastq_files
bamDir=$workDir/Bams
refGenome=/home/mjohnson/Desktop/Hemp_diversity_GS/Hemp_variant_calling/GCF_900626175.1_cs10_genomic.fna

# num of processors
PROCS=8

if [ ! -e $bamDir ] ; then mkdir $bamDir ; fi 


# Index reference genome for alignment
bwa index -p $dataDir/HempIndex $refGenome
samtools faidx $refGenome


# # pre-processing steps and initial variant calling
for sample in $(ls $dataDir | grep "_R2" | awk -F'_' '{print $1}'); do 
    
    # align
    bwa mem -t $PROCS $dataDir/HempIndex $dataDir/${sample}_R1_clipped_passed-re-filter.fastq.gz $dataDir/${sample}_R2_clipped_passed-re-filter.fastq.gz > $bamDir/$sample.aligned.sam
    
    # convert from sam to bam
    samtools view -b $bamDir/$sample.aligned.sam | samtools sort > $bamDir/$sample.aligned.bam
    
done
