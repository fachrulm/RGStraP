#!/bin/sh

## The script for PICARD's AddOrReplaceReadGroups. Please run 'addRG.sh' to apply the read group information to your samples, as this script applies to one sample.
## Please modify the path to your PICARD directory.

i=$1 #BAM file path
nodir=${i##*/} #get BAM file name
id=$2 #Read-Group ID
pl=$3 #Read-Group platform
sm=$4 #Read-Group sample name
lb=$5 #Read-Group library
pu=$6 #Read-Group platform unit

mkdir -p ./addRG/

picard AddOrReplaceReadGroups I=$i O=./addRG/${nodir%%_Aligned*}_withRG.bam RGID=$id RGPL=$pl RGSM=$sm RGLB=$lb RGPU=$pu

