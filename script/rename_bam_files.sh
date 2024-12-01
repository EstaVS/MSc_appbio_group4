#!/bin/bash

# Define the base directory
BASE_DIR="/scratch_tmp/grp/msc_appbio/group4_tmp/sorted_and_indexed_output"

# Loop through ERR4553381 to ERR4553384 directory
for SEQ_NUM in ERR4553381 ERR4553382 ERR4553383 ERR4553384; do
	DIR="${BASE_DIR}/${SEQ_NUM}"

    if [ -d "$DIR" ]; then
        echo "Processing directory: $DIR"
        
        # Define the old and new names for the BAM file
        OLD_BAM="${DIR}/${SEQ_NUM}_sorted.bam"
        NEW_BAM="${DIR}/${SEQ_NUM}_Aligned.sortedByCoord.out_sorted.bam"
        
        # Define the old and new names for the BAI file
        OLD_BAI="${DIR}/${SEQ_NUM}_sorted.bam.bai"
        NEW_BAI="${DIR}/${SEQ_NUM}_Aligned.sortedByCoord.out_sorted.bam.bai"
        
        # Rename the BAM file
        if [ -f "$OLD_BAM" ]; then
            mv "$OLD_BAM" "$NEW_BAM"
            echo "Renamed: $OLD_BAM -> $NEW_BAM"
        else
            echo "File not found: $OLD_BAM"
        fi
        
        # Rename the BAI file
        if [ -f "$OLD_BAI" ]; then
            mv "$OLD_BAI" "$NEW_BAI"
            echo "Renamed: $OLD_BAI -> $NEW_BAI"
        else
            echo "File not found: $OLD_BAI"
        fi
    else
        echo "Directory not found: $DIR"
    fi
done

echo "Renaming complete."
