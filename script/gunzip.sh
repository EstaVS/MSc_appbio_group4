#!/bin/bash

# Define source and target directories
SOURCE_DIR="../../raw_data_rna/fastq/"       # Directory containing the .gz files
TARGET_DIR="../../unzipped_data"        # Directory to save unzipped files

# Create the target directory if it doesn't exist
mkdir -p $TARGET_DIR

# Loop through all .gz files in the source directory
for gz_file in $SOURCE_DIR/*.gz; do
    # Extract the base filename without extension
    base_name=$(basename "$gz_file" .gz)

    # Unzip the file into the target directory
    echo "Unzipping $gz_file to $TARGET_DIR/$base_name"
    gunzip -c "$gz_file" > "$TARGET_DIR/$base_name"
done

echo "All files have been unzipped to $TARGET_DIR"
