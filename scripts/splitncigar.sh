#!/bin/bash

##JAVA VERSION 1.8 NEEDED, PLS SOURCE PROFILE

mkdir -p ./NCIGAR
bam=$1
file=${bam##*/}
genome=$2

#SplitNCigar
gatk --java-options "-Xmx20g -XX:ParallelGCThreads=12" SplitNCigarReads \
	-R $genome \
	-I $bam \
	-O ./NCIGAR/${file%%.bam}_splitN.bam

echo "Finished $file"


