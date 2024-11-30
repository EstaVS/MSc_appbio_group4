#!/bin/bash

#SBATCH --job-name=star_alignment     # Job name
#SBATCH --output=star_alignment.log   # Output log file
#SBATCH --error=star_alignment.err    # Error log file
#SBATCH --time=06:00:00               # Time limit (hh:mm:ss)
#SBATCH --cpus-per-task=4             # Number of CPU cores
#SBATCH --mem=32G                     # Memory pool for the job
#SBATCH --partition=msc_appbio          # Partition to submit to (adjust based on your cluster)
#SBATCH --mail-type=END,FAIL          # Notifications for job done & fail
#SBATCH --mail-user=k24023654@kcl.ac.uk  # Email for notification

# Parameters - customize these paths
INDEX_PATH="../../index/"  # Path to STAR indexed genome
REF_GENOME_DIR="/scratch_tmp/grp/msc_appbio/group4_tmp/ref_genome/GCA_000269885.1_ASM26988v1_genomic.fna"    # Path to reference genome
FASTQ_DIR="../../raw_data_rna/fastq/"
OUTPUT_DIR="../../star_output"                     # Directory for STAR output
THREADS=4                                    # Number of threads for STAR

# load star
module load star/2.7.10b-gcc-13.2.0

# Create output directory if it doesn't exist
mkdir -p $OUTPUT_DIR

# Loop through each fastq.gz file in the raw_data_rna/fastq directory
for FASTQ in $FASTQ_DIR/*.fastq.gz;do
	BASENAME=$(basename $FASTQ .fastq.gz)

	# Run STAR
	echo "Running STAR alignment for $BASENAME..."

	STAR --genomeDir $INDEX_PATH \
    	--readFilesIn $FASTQ \
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
   	 --outFileNamePrefix ${OUTPUT_DIR}/${BASENAME}_ \
   	 --outSAMattributes NH HI AS NM MD \
   	 --outSAMstrandField intronMotif \
   	 --outSAMtype BAM SortedByCoordinate \
	--quantMode TranscriptomeSAM GeneCounts \
    	--sjdbScore 1 \
    	--limitBAMsortRAM 8000000000

	echo "Alignment completed for $BASENAME. Output saved to ${OUTPUT_DIR}/${BASENAME}_"
done
