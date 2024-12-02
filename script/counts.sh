#!/bin/bash 

#Define paths
baseDir="/scratch_tmp/grp/msc_appbio/group4_tmp"
bamDir="${baseDir}/sorted_and_indexed_output"
gtf_file="${baseDir}/ref_genome/genomic.gtf"
outputDir="${baseDir}/counts"

#Create output directory if it doesn't exist 
mkdir -p "$outputDir"

#Loop through all BAM files in the subdirectories matching ERR45533*
for bamfile in "$bamDir"/ERR45533*/*sortedByCoord.out_sorted.bam; do

	#Extract sample name from, the BAM file.
	sampleName=$(basename "$bamfile" .bam)

	#Define the output counts file 
	outputCounts="${outputDir}/counts_${sampleName}.txt"

	#Run HTSeq-Count for the BAM files, annotating it with GTF file
	echo "Running HTSeq-counts for sample" "$sampleName"
	htseq-count -f bam -r pos -s no "$bamfile" "$gtf_file" > "$outputCounts"

	 # Check if HTSeq-Count ran successfully
        if [[ $? -eq 0 ]]; then
            echo "Counts for $sampleName directed to  $outputCounts"
        else
            echo "Error: HTSeq-Count failed for $sampleName. Check input files."
        fi
done


echo "Gene counting with HTSeq is completed for all samples. Results are in $outputDir."
