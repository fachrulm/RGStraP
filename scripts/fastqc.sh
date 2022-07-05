#!/bin/bash
# SLURM_dl.sh

p1=$1

mkdir -p ./fastqc
mkdir -p ./logs

fastqc $p1 -f fastq -o ./fastqc
fastqc ${p1%1.fastq.gz}2.fastq.gz -f fastq -o ./fastqc

