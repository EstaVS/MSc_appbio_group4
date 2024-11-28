#!/bin/bash

# Define source and target directories
SOURCE_DIR="../../processed_data"       # Directory containing the .zip files
TARGET_DIR="../../unzipped_data"        # Directory to save unzipped contents

# Create the target directory if it doesn't exist
mkdir -p $TARGET_DIR

# Loop through all .zip files in the source directory
for zip_file in $SOURCE_DIR/*.zip; do
    # Extract the base filename without extension
    base_name=$(basename "$zip_file" .zip)

    # Create a subdirectory for each unzipped file
    output_subdir="$TARGET_DIR/$base_name"
    mkdir -p "$output_subdir"

    # Unzip the file into the target directory
    echo "Unzipping $zip_file to $output_subdir"
    unzip -q "$zip_file" -d "$output_subdir"
done

echo "All FastQC reports have been unzipped to $TARGET_DIR"
