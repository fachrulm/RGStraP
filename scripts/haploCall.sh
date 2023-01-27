#!/bin/bash

##JAVA VERSION 1.8 NEEDED, PLS SOURCE PROFILE

mkdir -p ./vars

bam=$1
file=${bam##*/}
genome=$2

#HaplotypeCaller
gatk HaplotypeCaller \
	-R $genome \
	-I $bam \
	-ERC GVCF \
	-O ./vars/${file%%_final*}_first.g.vcf.gz


