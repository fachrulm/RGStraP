#!/bin/sh

## The script for PICARD's AddOrReplaceReadGroups. Please run 'addRG.sh' to apply the read group information to your samples, as this script applies to one sample.
## Please modify the path to your PICARD directory.

i=$1 #BAM file path
nodir=${i##*/} #get BAM file name
PICARD=$2
id=$3 #Read-Group ID
pl=$4 #Read-Group platform
sm=$5 #Read-Group sample name
lb=$6 #Read-Group library
pu=$7 #Read-Group platform unit

mkdir -p ./addRG/

java -jar $PICARD AddOrReplaceReadGroups I=$i O=./addRG/${nodir%%_Aligned*}_withRG.bam RGID=$id RGPL=$pl RGSM=$sm RGLB=$lb RGPU=$pu

