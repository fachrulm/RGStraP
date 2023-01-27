#!/bin/bash

##JAVA VERSION 1.8 NEEDED, PLS SOURCE PROFILE

mkdir -p ./VCF

list=$1
genome=$2
int=$3

#Merge GVCFs
gatk CombineGVCFs \
	-R $genome \
	$(cat $list) \
	-O ./VCF/Merged.g.vcf.gz

#Genotype GVCFs
gatk GenotypeGVCFs \
	-R $genome \
	--intervals $int \
	-V ./VCF/Merged.g.vcf.gz \
	-O ./VCF/Genotyped_raw.vcf.gz
