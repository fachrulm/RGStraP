#!/bin/bash

##JAVA VERSION 1.8 NEEDED, PLS SOURCE PROFILE
#source ~/.bash_profile

mkdir -p ./VCF

vcf=$1
genome=$2
gatk4=$3

#Genotype GVCFs
$gatk4 VariantFiltration \
	-R $genome \
	-V $vcf \
	--filter-expression "QD < 2.0" --filter-name "QDfilterLessThan2" \
	--filter-expression "MQ < 40.0" --filter-name "MQfilterLessThan40" \
	--filter-expression "FS > 60.0" --filter-name "FSfilterMoreThan60" \
	--filter-expression "MQRankSum < -12.5" --filter-name "MQRankSumfilterLessThan-12.5" \
	--filter-expression "ReadPosRankSum < -8.0" --filter-name "ReadPosfilterLessThan-8.0" \
	--genotype-filter-expression "DP < 2" --genotype-filter-name "DP3" \
       	--set-filtered-genotype-to-no-call true \
	-O ${vcf%%_raw*}_filterKeep.vcf.gz

