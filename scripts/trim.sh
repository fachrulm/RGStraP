#!/bin/bash

p1=$1

mkdir -p ./trimmed

echo "Trimming ${p1##*/}"
trim_galore -q 20 --length 20 -o ./trimmed/ --paired $p1 ${p1%1.fastq.gz}2.fastq.gz

