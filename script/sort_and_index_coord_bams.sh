#!/bin/bash
#SBATCH --job-name=sort_index_bam  # Job name
#SBATCH --output=sort_index_%j.log # Standard output log file
#SBATCH --error=sort_index_%j.err  # Standard error log file
#SBATCH --time=02:00:00            # Time limit (hh:mm:ss)
#SBATCH --mem=8G                   # Memory limit
#SBATCH --ntasks=1                 # Number of tasks
#SBATCH --cpus-per-task=1          # Number of CPU cores
#SBATCH --partition=msc_appbio     # Partition to use

# Load samtools module
module load samtools/1.17-gcc-13.2.0-python-3.11.6

# Base directory containing the subfolders with BAM files
BASE_DIR="/scratch_tmp/grp/msc_appbio/group4_tmp/sorted_and_indexed_output"

echo "Starting sorting and indexing of sortedByCoord BAM files..."

# Loop through each folder and process the BAM file
for SEQ_DIR in "$BASE_DIR"/*; do
    # Check if SEQ_DIR is a directory
    if [ -d "$SEQ_DIR" ]; then
        # Find the BAM file in the directory
        BAM_FILE=$(find "$SEQ_DIR" -maxdepth 1 -type f -name "*_Aligned.sortedByCoord.out.bam")
        
        # Skip if no BAM file is found
        if [ -z "$BAM_FILE" ]; then
            echo "No BAM file found in $SEQ_DIR, skipping..."
            continue
        fi

        # Extract file and output names
        BAM_FILENAME=$(basename "$BAM_FILE")       # Original BAM file name
        OUTPUT_BAM="$SEQ_DIR/$BAM_FILENAME"       # Sorted BAM file

        echo "Processing file: $BAM_FILE"

        # Sort BAM file by coordinate (in case further sorting is needed)
        samtools sort -o "$OUTPUT_BAM" "$BAM_FILE"

        # Index the sorted BAM file
        samtools index "$OUTPUT_BAM"

        echo "Sorted and indexed: $OUTPUT_BAM"
    fi
done

echo "All BAM files have been processed."

