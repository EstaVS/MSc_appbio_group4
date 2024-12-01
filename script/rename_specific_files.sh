#!/bin/bash

# Define the base directory
BASE_DIR="/scratch_tmp/grp/msc_appbio/group4_tmp/sorted_and_indexed_output"

# Define the range of directories to process
for i in {381..384}; do
    DIR="$BASE_DIR/ERR4553$i"
    
    # Check if the directory exists
    if [ -d "$DIR" ]; then
        echo "Processing directory: $DIR"
        cd "$DIR" || continue
        
        # Rename files by removing ".fastq.gz"
        for file in *.fastq.gz*; do
            if [ -f "$file" ]; then
                new_name="${file/.fastq.gz/}"
                mv "$file" "$new_name"
                echo "Renamed: $file -> $new_name"
            fi
        done
    else
        echo "Directory not found: $DIR"
    fi
done

echo "Renaming complete for specified directories."
