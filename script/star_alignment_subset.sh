#!/bin/bash

#SBATCH --job-name=star_alignment_subset     # Job name
#SBATCH --output=star_alignment_subset.log   # Output log file
#SBATCH --error=star_alignment_subset.err    # Error log file
#SBATCH --time=06:00:00               # Time limit (hh:mm:ss)
#SBATCH --cpus-per-task=4             # Number of CPU cores
#SBATCH --mem=32G                     # Memory pool for the job
#SBATCH --partition=msc_appbio          # Partition to submit to (adjust based on your cluster)
#SBATCH --mail-type=END,FAIL          # Notifications for job done & fail
#SBATCH --mail-user=k24023654@kcl.ac.uk  # Email for notification

# star_alignment_subset.sh
module load star/2.7.10b-gcc-13.2.0 

# Directory paths
RAW_DATA_DIR="/scratch_tmp/grp/msc_appbio/group4_tmp/raw_data_rna/fastq"
OUTPUT_DIR="/scratch_tmp/grp/msc_appbio/group4_tmp/star_output_subset"
INDEX_PATH="../../index/"  # Path to STAR indexed genome 
THREADS=4

# Specify the sample IDs to process
SAMPLES=("ERR4553381.fastq.gz" "ERR4553382.fastq.gz" "ERR4553383.fastq.gz" "ERR4553384.fastq.gz")

# Iterate over each sample ID
for SAMPLE in "${SAMPLES[@]}"; do

    #Full path
    FULL_PATH="${RAW_DATA_DIR}/${SAMPLE}"

    #Extract the basename without extension
    SEQ_NUM=$(basename "$SAMPLE" .fastq.gz)
    echo "Processing $SEQ_NUM..."

    STAR  --genomeDir $INDEX_PATH \
        --readFilesIn $FULL_PATH \
        --readFilesCommand zcat \
         --runThreadN $THREADS \
         --genomeLoad NoSharedMemory \
         --outFilterMultimapNmax 20 \
         --alignSJoverhangMin 8 \
         --alignSJDBoverhangMin 1 \
         --outFilterMismatchNmax 999 \
         --outFilterMismatchNoverReadLmax 0.04 \
         --alignIntronMin 20 \
         --alignIntronMax 1000000 \
         --alignMatesGapMax 1000000 \
         --outSAMheaderHD @HD VN:1.4 SO:coordinate \
         --outSAMunmapped Within \
         --outFilterType BySJout \
         --outFileNamePrefix ${OUTPUT_DIR}/${SAMPLE}_ \
         --outSAMattributes NH HI AS NM MD \
         --outSAMstrandField intronMotif \
         --outSAMtype BAM SortedByCoordinate \
        --quantMode TranscriptomeSAM GeneCounts \
        --sjdbScore 1 \
        --limitBAMsortRAM 8000000000

        echo "Alignment completed for $BASENAME. Output saved to ${OUTPUT_DIR}/${SEQ_NUM}_"
done

