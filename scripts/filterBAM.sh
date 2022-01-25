#!/bin/bash

## Filter reads with low MAPQ (<20) with samtools. 


bam=$1 #BAM file
sname=${bam##*/} #without directory
samtools=$2

$samtools view -bq 20 $bam > ${bam%.*}_filtered.bam

