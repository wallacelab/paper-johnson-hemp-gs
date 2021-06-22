# Analysis pipeline for Johnson et al. 2021

This is the analysis pipeline for Johnson et al 2021, "Genomic and Chemical Diversity of Commercially Available Industrial Hemp Accessions." It contains all analysis scripts starting from turning aligned BAM files into genotypes. It also includes key intermediate files so that users can pick up partway instead of having to run the entire thing from scratch.

# Setting up the Conda environment
This pipieline uses a Conda environment for reproducibility.  To set up an identical environment:
1. Install [Anaconda](https://docs.anaconda.com/anaconda/install/) (I recommend the Miniconda version)
2. Follow the instructions in 0_CreateCondaEnvironment.sh to initialize an identical environment

# Running the analysis
After setting up your conda enviroment (above), you can run 0_RerunHempGenotyping.sh. The script assumes you have already aligned the FASTQ files to a genome and gotten BAM files out. If you don't want to rerun the entire pipeline from scratch, you can start any of the key intermediate files included in this repo.
* FASTQ files from the paper (adaptor-trimmed and restriction-fragment-verified): [NCBI Bioproject PRJNA707556](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA707556)
* Hemp genome: [CBDRx version GCF_900626175.1](https://www.ncbi.nlm.nih.gov/assembly/GCF_900626175.1/)

# Key files in this repository
* Genotypes: 1g_hemp_filtered.sorted.vcf.gz
* Keyfile linking sample number to accession name: 0_HDGS_keyfile.txt
* HPLC results from each sample: 0_HDGS_HPLC_Results.csv
