#!/bin/sh

## Add read group information to BAM files, which are necessary during GATK's variant calling process. Fastq files do not natively have this information, hence we are required to manually curate the information.
## Running this requires a list of BAM files as well as a space-delimited list (with no header) containing read group information for each sample, which include (in order:
## RGID = Read-Group ID
## RGPL = Read-Group platform
## RGSM = Read-Group sample name
## RGLB = Read-Group library
## RGPU = Read-Group platform unit
## For more information: https://gatk.broadinstitute.org/hc/en-us/articles/360037226472-AddOrReplaceReadGroups-Picard-

bam=$1 #BAM file
sname=${bam##*/} #without directory

meta=$2 #read group information, file needs to be ordered the same with list of BAM files
RGinfo=$(grep "${sname%_Al*}" $meta | sed 's/.*.bam //g') #get relevant RG info per sample

addRG=$3 #script to run the operation per sample

mkdir -p ./addRG/

echo "Adding read group information to ${sname%%_Aligned*}"
sh $addRG $bam $RGinfo
echo "successful!"

