#!/bin/bash

##JAVA VERSION 1.8 NEEDED, PLS SOURCE PROFILE
#source ~/.bash_profile

mkdir -p ./VCF

list=$1
genome=$2
int=$3

#Merge GVCFs
gatk4 CombineGVCFs \
	-R $genome \
	$(cat $list) \
	-O ./VCF/Merged.g.vcf.gz

#Genotype GVCFs
gatk4 GenotypeGVCFs \
	-R $genome \
	--intervals $int \
	-V ./VCF/Merged.g.vcf.gz \
	-O ./VCF/Genotyped_raw.vcf.gz
