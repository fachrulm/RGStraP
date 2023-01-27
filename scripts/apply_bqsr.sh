#!/bin/bash

##JAVA VERSION 1.8 NEEDED, PLS SOURCE PROFILE

mkdir -p ./bqsr
bam=$1
file=${bam##*/}
genome=$2
rectab=$3

#Apply BQSR
gatk ApplyBQSR \
	-I $bam \
	-R $genome \
	--bqsr-recal-file $rectab \
	-O ./bqsr/${file%%_withRG*}_final.bam

