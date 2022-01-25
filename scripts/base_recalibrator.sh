#!/bin/bash

##JAVA VERSION 1.8 NEEDED, PLS SOURCE PROFILE
#source ~/.bash_profile

mkdir -p ./bqsr
bam=$1
file=${bam##*/}
genome=$2
gatk4=$3
indel1=$4
indel2=$5

#BaseRecalibrator
$gatk4 BaseRecalibrator \
	-I $bam \
	-R $genome \
	--known-sites $indel1 \
	--known-sites $indel2 \
	-O ./bqsr/${file%%_withRG*}_recal.table

