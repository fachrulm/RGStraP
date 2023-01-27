#!/bin/bash

##JAVA VERSION 1.8 NEEDED, PLS SOURCE PROFILE

mkdir -p ./VCF

vcf=$1
genome=$2

#Genotype GVCFs
gatk SelectVariants \
	-R $genome \
	-V $vcf \
	--exclude-filtered --exclude-non-variants \
	-O ${vcf%%_filterKeep*}_filterOut.vcf.gz
