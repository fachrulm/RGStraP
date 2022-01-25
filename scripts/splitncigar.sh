#!/bin/bash

##JAVA VERSION 1.8 NEEDED, PLS SOURCE PROFILE
#source ~/.bash_profile

mkdir -p ./NCIGAR
bam=$1
file=${bam##*/}
genome=$2
gatk4=$3

#SplitNCigar
$gatk4 --java-options "-Xmx20g -XX:ParallelGCThreads=12" SplitNCigarReads \
	-R $genome \
	-I $bam \
	-O ./NCIGAR/${file%%.bam}_splitN.bam

echo "Finished $file"


