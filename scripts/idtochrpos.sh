#!/bin/bash

bim=$1

cp $bim ${bim}_bck
awk '{if($2 == ".") {print $1,"\t",$1":"$4,"\t",$3,"\t",$4,"\t",$5,"\t",$6} else {print $0}}' ${1}_bck > $1
