#!/bin/bash

# Parameters - customize these paths
INDEX_PATH="../../index/SAindex"
REF_GENOME_DIR="/scratch_tmp/grp/msc_appbio/group4_tmp/ref_genome/GCA_000269885.1_ASM26988v1_genomic.fna"    # Path to STAR indexed genome
FASTQ="../../raw_data_rna/fastq/ERR4553381.fastq.gz"
OUTPUT_DIR="../../star_output_test381"                     # Directory for STAR output
THREADS=4                                    # Number of threads for STAR

# load star
module load star/2.7.10b-gcc-13.2.0

# Create output directory if it doesn't exist
mkdir -p $OUTPUT_DIR

# Run STAR
echo "Running STAR alignment..."

STAR --genomeDir $INDEX_PATH \
    --readFilesIn $FASTQ \
    --readFilesCommand zcat --runThreadN 4 \
    --genomeLoad NoSharedMemory --outFilterMultimapNmax 20 \
    --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 \
    --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04 \
    --alignIntronMin 20 --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 --outSAMheaderHD @HD VN:1.4 SO:coordinate \
    --outSAMunmapped Within --outFilterType BySJout --outFileNamePrefix ${OUTPUT_DIR} \
    --outSAMattributes NH HI AS NM MD --outSAMstrandField intronMotif \
    --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM GeneCounts \
    --sjdbScore 1 --limitBAMsortRAM 8000000000 \

echo "Alignment completed. Output saved to ${OUTPUT_DIR}/"
