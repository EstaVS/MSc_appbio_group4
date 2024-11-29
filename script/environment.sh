#!/bin/bash

echo "start of the pipeline"

# Add conda channels
conda config --add channels bioconda
conda config --add channels conda-forge

# Create a new conda environment with specific versions
conda create -n bioenv cutadapt=2.3 faqcs=2.08 star= 2.5.2a fdrtool= 1.2.15 htseq 0.9.1 -y

#Load Fastqc

module load fastqc

# Print a message about how to activate the environment
echo "To activate your environment, run: conda activate bioenv" 


echo "end of the pipeline"
