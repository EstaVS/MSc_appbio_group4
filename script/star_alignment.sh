#!/bin/bash

# Parameters - customize these paths
REF_GENOME_DIR="/scratch_tmp/grp/msc_appbio/group4_tmp/ref_genome/GCA_000269885.1_ASM26988v1_genomic.fna"    # Path to STAR indexed genome
FASTQ="../../raw_data_rna/fastq/ERR4553381.fastq.gz"
OUTPUT_DIR="../../star_output"                     # Directory for STAR output
THREADS=4                                    # Number of threads for STAR

# load star
module load star/2.7.10b-gcc-13.2.0

# Create output directory if it doesn't exist
mkdir -p $OUTPUT_DIR

# Run STAR
echo "Running STAR alignment..."
STAR \
        --genomeDir $REF_GENOME_DIR \
        --readFilesIn $FASTQ \
        --readFilesCommand zcat \
        --runThreadN $THREADS \
        --outFileNamePrefix ${OUTPUT_DIR}/ \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMattributes Standard \
        --quantMode TranscriptomeSAM GeneCounts\
        --peOverlapNbasesMin 10

echo "Alignment completed. Output saved to ${OUTPUT_DIR}/"
