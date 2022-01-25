#!/bin/bash

file=$1   #prefix name of the plink (.bim/.bed/.fam) file
output=$2 #name of output file
#value 0 means that the filter is not used
maf=$3  #ex. 0.01
hwe=$4  #ex. 0.000001
geno=$5 #ex. 0.1  (missingnes per-variant)
mind=$6 #ex. 0.1  (missigness per-sample)
plink1_9=$7 #path to plink1.9


tmp_file="$file""_tmp"
$plink1_9 --bfile $file --maf $maf --hwe $hwe --geno $geno --allow-no-sex --make-bed --out $tmp_file
$plink1_9 --bfile $tmp_file --mind $mind --allow-no-sex --make-bed --out $output
rm $tmp_file.*


