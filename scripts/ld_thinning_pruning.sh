#!/bin/bash

file1=$1
output=$2
win=$3
step=$4
r2=$5
excl=$6
plink1_9=$7

$plink1_9 --bfile $file1 --indep-pairwise $win $step $r2 --exclude range ${excl} --allow-no-sex
$plink1_9 --bfile $file1 --extract plink.prune.in --make-bed --allow-no-sex --out $output 
rm plink.*
