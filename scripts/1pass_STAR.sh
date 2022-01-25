#!/bin/bash
## Perform first-pass mapping using STAR.
## If you encounter the 'not enough memory for BAM sorting' error, add the --limitBAMsortRAM line with the recommended amount.

p1=$1
STAR_genome=$2 #path to STAR-indexed genome build
mkdir -p ./1_mapped

sname=${p1##*/}
echo "Mapping ${sname%_1_val_1.fq.gz}"
STAR --runThreadN 12 --runMode alignReads \
	--readFilesCommand zcat \
	--outSAMtype BAM SortedByCoordinate \
	--outFileNamePrefix ./1_mapped/${sname%_1_val_1_clumped.fq.gz}_ \
	--outSAMmapqUnique 60 \
	--genomeDir ${STAR_genome} \
	--readFilesIn $i ${p1%_1_val_1_clumped.fq.gz}_2_val_2_clumped.fq.gz
