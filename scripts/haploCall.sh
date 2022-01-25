#!/bin/bash

##JAVA VERSION 1.8 NEEDED, PLS SOURCE PROFILE
#source ~/.bash_profile

mkdir -p ./vars

bam=$1
file=${bam##*/}
genome=$2
gatk4=$3

#HaplotypeCaller
$gatk4 HaplotypeCaller \
	-R $genome \
	-I $bam \
	-ERC GVCF \
	-O ./vars/${file%%_final*}_first.g.vcf.gz


