#! /usr/bin/env bash 


#################
# OVERVIEW
#################

# To make this work reproducible, most of the analyses were done using Conda environments. 
# This file contains the code necessary to recreate those environments, along with how
# they were created in the first place. (Mostly for documentations sake)  

# Environment name - This is what the conda environment will be called. You can change this if you wish.
conda_env=proj-hemp-genotyping-2020



#################
# RECREATE THE ENVIRONMENT USED IN THE ANALYSIS
#################

# To recreate the environment used in this analysis, first make sure you have some version of anaconda installed. 
# Then run the below command. This will use the included conda_environment.txt specifications to exactly recreate the environment

conda create --name $conda_env --file conda_environment.txt



#################
# LOADING THE ENVIRONMENT
#################

# Once you have created the environment, you can activate (load) it by running the following command.
# If you're doing this on the command line, make sure you type out the full name of the environment, 
# since the $conda_env variable won't be active in a new shell.
# conda activate $conda_env


# # If you want to run the above command inside a script (like I prefer), you need to include the following line of code first:
# . $(conda info --root)/etc/profile.d/conda.sh 


# So in a script, you'd probably want to include these three lines of code near the begining (be sure to uncomment them:
# conda_env=proj-hemp-genotyping-2020
# . $(conda info --root)/etc/profile.d/conda.sh 
# conda activate $conda_env



#################
# ORIGINAL COMMANDS - DO NOT USE
#################

# These are the commands originally used to create the above environment (and export it to conda_environment.txt).
# You should not need to run these commands, but they are included to document what was done.

# # Create environment; scipy added later, hence in another command
# conda create -n $conda_env bcftools=1.9 vcftools=0.1.16 tassel=5.2.40 r-base=3.6.2 r-essentials=3.6
# conda activate $conda_env
# conda install scipy


# # Export environment for others to load
# conda list --explicit > conda_environment.txt









