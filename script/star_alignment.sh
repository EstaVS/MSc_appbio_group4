#!/bin/bash

# Parameters - customize these paths
REF_GENOME_DIR="path_to_reference_genome"    # Path to STAR indexed genome
FASTQ_1="reads_1.fastq.gz"                   # Input FASTQ file (R1)
FASTQ_2="reads_2.fastq.gz"                   # Input FASTQ file (R2) for paired-end; omit for single-end
OUTPUT_DIR="../star_output"                     # Directory for STAR output
THREADS=4                                    # Number of threads for STAR

# Create output directory if it doesn't exist
mkdir -p $OUTPUT_DIR

# Run STAR
echo "Running STAR alignment..."
if [[ -f "$FASTQ_2" ]]; then
    STAR \
        --genomeDir $REF_GENOME_DIR \
        --readFilesIn $FASTQ_1 $FASTQ_2 \
        --readFilesCommand zcat \
        --runThreadN $THREADS \
        --outFileNamePrefix ${OUTPUT_DIR}/ \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMattributes Standard \
        --quantMode TranscriptomeSAM GeneCounts
else
    STAR \
        --genomeDir $REF_GENOME_DIR \
        --readFilesIn $FASTQ_1 \
        --readFilesCommand zcat \
        --runThreadN $THREADS \
        --outFileNamePrefix ${OUTPUT_DIR}/ \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMattributes Standard \
        --quantMode TranscriptomeSAM GeneCounts
fi

echo "Alignment completed. Output saved to ${OUTPUT_DIR}/"
