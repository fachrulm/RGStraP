#!/bin/bash

##JAVA VERSION 1.8 NEEDED, PLS SOURCE PROFILE
#source ~/.bash_profile

mkdir -p ./VCF

vcf=$1
genome=$2
gatk4=$3

#Genotype GVCFs
$gatk4 SelectVariants \
	-R $genome \
	-V $vcf \
	--exclude-filtered --exclude-non-variants \
	-O ${vcf%%_filterKeep*}_filterOut.vcf.gz
