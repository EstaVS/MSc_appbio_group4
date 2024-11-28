#!/bin/bash

# Input files
FASTQC_ZIP="../processed_data/ERR4553381_fastqc.zip"  # Replace with your FastQC .zip file
STAR_LOG="Log.final.out"         # Replace with your STAR alignment log file

# Function to extract FastQC summary
check_fastqc() {
    echo "Checking FastQC results..."
    unzip -p $FASTQC_ZIP "*/summary.txt" | awk '
    BEGIN {print "FastQC Summary:"}
    {print $2, "-", $1}
    '
}

# Function to analyze STAR alignment log
check_star_alignment() {
    echo "Checking STAR alignment results..."
    awk '
    /Number of input reads/ {print "Input reads:", $NF}
    /Uniquely mapped reads %/ {print "Uniquely mapped reads (%):", $NF}
    /% of reads mapped to multiple loci/ {print "Mapped to multiple loci (%):", $NF}
    /% of reads mapped to too many loci/ {print "Mapped to too many loci (%):", $NF}
    /% of reads unmapped: too many mismatches/ {print "Unmapped due to mismatches (%):", $NF}
    /% of reads unmapped: too short/ {print "Unmapped due to short reads (%):", $NF}
    ' $STAR_LOG
}

# Run checks
if [[ -f $FASTQC_ZIP && -f $STAR_LOG ]]; then
    check_fastqc
    echo ""
    check_star_alignment
else
    echo "Error: Missing input files. Ensure $FASTQC_ZIP and $STAR_LOG exist."
fi
