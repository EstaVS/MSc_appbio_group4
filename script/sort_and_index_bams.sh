#!/bin/bash
#SBATCH --job-name=sort_index_bam   # Job name
#SBATCH --output=sort_index_bam.out # Output log file
#SBATCH --error=sort_index_bam.err  # Error log file
#SBATCH --ntasks=1                  # Number of tasks (1 job)
#SBATCH --cpus-per-task=4           # Number of CPUs per task
#SBATCH --mem=16G                   # Memory per node
#SBATCH --time=2:00:00              # Time limit (hh:mm:ss)
#SBATCH --partition=msc_appbio      # Partition/queue name
#SBATCH --mail-type=END,FAIL        # Notifications for job done & fail

# Load necessary modules
module load samtools/1.17-gcc-13.2.0-python-3.11.6

# Input and output directories
#INPUT_DIR="../../star_output"        # Directory containing BAM files
INPUT_DIR="../../star_output_subset"  # Directory containing BAM files realigned
#OUTPUT_DIR="../../sorted_and_indexed_output"  # Directory to store organized files
OUTPUT_DIR="../../sorted_and_indexed_output_subset"  #Directory to store organized realigned files
THREADS=4                      # Number of threads for sorting

# Create output directory if it does not exist
mkdir -p $OUTPUT_DIR

# Process each BAM file
for BAM_FILE in ${INPUT_DIR}/*Aligned.sortedByCoord.out.bam; do
    # Extract sequence number (e.g., ERR4553381)
    SEQ_NUM=$(basename $BAM_FILE | cut -d_ -f1)

    # Create a directory for this sequence number
    SEQ_DIR=${OUTPUT_DIR}/${SEQ_NUM}
    mkdir -p $SEQ_DIR

    # Output file paths
    SORTED_BAM=${SEQ_DIR}/${SEQ_NUM}_sorted.bam
    BAM_INDEX=${SEQ_DIR}/${SEQ_NUM}_sorted.bam.bai

    echo "Processing $SEQ_NUM..."

    # Sort BAM file
    samtools sort -@ $THREADS -o $SORTED_BAM $BAM_FILE

    # Index the sorted BAM file
    samtools index $SORTED_BAM

    # Move original and log files into the sequence directory
    mv ${INPUT_DIR}/${SEQ_NUM}* ${SEQ_DIR}/

    echo "Finished processing $SEQ_NUM. Outputs saved to $SEQ_DIR"
done

echo "All BAM files processed and organized."
