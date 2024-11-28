#!/bin/bash


echo "start of the pipeline"

#Load fastqc

module load fastqc 

#Define base directory and the results directory 

baseDirectory="/scratch_tmp/grp/msc_appbio/group4_tmp/raw_data_rna/fastq"
resultsDirectory="/scratch_tmp/grp/msc_appbio/group4_tmp/processed_data"

#Create the results directory if it does not exist

mkdir -p "$resultsDirectory"

#Run FASTQC on all .fastq.gz files in the base directory 

fastqc -o "$resultsDirectory" -t 4 "$baseDirectory"/*.fastq.gz

echo "end of the pipeline"
