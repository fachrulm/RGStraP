#!/bin/bash

p1=$1
ext=${p1##*1.}

mkdir -p ./trimmed

echo "Trimming ${p1##*/}"
#echo "extension ${ext}"
#echo "${p1%1.fastq*}2.${ext}"
trim_galore -q 20 --length 20 --gzip -o ./trimmed/ --paired $p1 ${p1%1.fastq*}2.${ext}

