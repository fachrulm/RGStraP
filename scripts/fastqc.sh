#!/bin/bash
# SLURM_dl.sh

p1=$1
#fastqc=$2

mkdir -p ./fastqc

fastqc $p1 -f fastq -o ./fastqc
fastqc ${p1%1.fastq.gz}2.fastq.gz -f fastq -o ./fastqc

