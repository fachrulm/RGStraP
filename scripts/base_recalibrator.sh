#!/bin/bash

##JAVA VERSION 1.8 NEEDED, PLS SOURCE PROFILE

mkdir -p ./bqsr
bam=$1
file=${bam##*/}
genome=$2
indel1=$3
indel2=$4

#BaseRecalibrator
gatk BaseRecalibrator \
	-I $bam \
	-R $genome \
	--known-sites $indel1 \
	--known-sites $indel2 \
	-O ./bqsr/${file%%_withRG*}_recal.table

