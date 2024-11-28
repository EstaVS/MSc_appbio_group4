#!/bin/bash 

#Loading STAR
module load star


STAR --runThreadN 8 --runMode genomeGenerate --genomeDir /scratch_tmp/grp/msc_appbio/group4_tmp/ref_genome --genomeFastaFiles /scratch_tmp/grp/msc_appbio/group4_tmp/ref_genome/GCA_000269885.1_ASM26988v1_genomic.fna --sjdbGTFfile /scratch_tmp/grp/msc_appbio/group4_tmp/ref_genome/genomic.gtf 
