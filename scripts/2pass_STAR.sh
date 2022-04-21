#!/bin/bash
## Perform second-pass mapping using STAR, taking into account splice junction information ("SJ.out.tab" files) from first-pass mapping.
## If you encounter the 'not enough memory for BAM sorting' error, add the --limitBAMsortRAM line with the recommended amount.

p1=$1
STAR_genome=$2 #path to STAR-indexed genome build
sjfile=$3 #splice junction file from first-pass mapping
mkdir -p ./2_mapped

sname=${p1##*/}
echo "Mapping ${sname%_1_val_1.fq.gz}"
STAR --runThreadN 12 --runMode alignReads \
	--readFilesCommand zcat \
	--outSAMtype BAM SortedByCoordinate \
	--outFileNamePrefix ./2_mapped/${sname%_1_val_1_clumped.fq.gz}_ \
	--genomeDir ${STAR_genome} \
	--outSAMmapqUnique 60 \
	--readFilesIn $i ${p1%_1_val_1_clumped.fq.gz}_2_val_2_clumped.fq.gz \
	--sjdbFileChrStartEnd ${sjfile}
