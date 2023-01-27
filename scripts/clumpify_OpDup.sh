#!/bin/bash
## Warning: please adjust the 'dupedist' value according to your sequencing platform.
## Recommendations from bbmap's author:
## Nextseq	40 (and spany=t)
## HiSeq 1T	40
## HiSeq 2500	40
## HiSeq 3k/4k	2500
## NovaSeq	12000

p1=$1

mkdir -p ./no_OpDup

sname=${p1##*/}
echo "Removing optical duplicate from ${sname%%_1_val_1.fq.gz} samples"
clumpify.sh -Xmx55G in1=$p1 in2=${p1%%_1_val_1.fq.gz}_2_val_2.fq.gz \
	out1=./no_OpDup/${sname%.fq.gz}_clumped.fq.gz out2=./no_OpDup/${sname%%_1_val_1.fq.gz}_2_val_2_clumped.fq.gz \
	dedupe allduplicates optical=t dupedist=40 subs=0 spantiles=f
