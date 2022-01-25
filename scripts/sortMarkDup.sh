#!/bin/bash
#source /usr/bin/java

# -I input.bam #the input file has read groups added and must be coordinate sorted
# -O output.bam file with duplicates marked
# --CREAT_INDEX set this to true to get ".bai" index file
# --VALIDATION_STRINGENCY
# -M #a file called metrics.txt containing duplicate read data metrics and a .bai index file

mkdir -p ./nodup

bam=$1
PICARD=$2

#Get file name without path
file=${bam##*/}

#SortSam to sort BAM files
java -jar $PICARD SortSam \
	-I $bam \
        -O ./nodup/temp_${file%%.bam}_marked.bam \
	--SORT_ORDER coordinate

#Picard MarkDuplicates to mark and remove duplicates
java -jar $PICARD MarkDuplicates \
	-I ./nodup/temp_${file%%.bam}_marked.bam \
	-O ./nodup/${file%%.bam}_marked.bam \
	--REMOVE_SEQUENCING_DUPLICATES true \
	--CREATE_INDEX true \
	--VALIDATION_STRINGENCY SILENT \
	-M ./nodup/${file%%.bam}_marked_dup_metrics.txt
	
#Remove temp file
rm ./nodup/temp_${file%%.bam}_marked.bam

