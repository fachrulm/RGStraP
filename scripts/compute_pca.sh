#!/bin/bash
#Partition: 'sysgen' (<24h) or 'sysgen_long' (>24h)
#SBATCH -p sysgen
#Give the job a name:
#SBATCH --job-name="pca"
# Maximum number of tasks/CPU cores used by the job:
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
# The amount of memory in megabytes per process in the job (try to put the smaller file first):
#SBATCH --mem=30000
# The maximum running time of the job in days-hours:mins:sec
#SBATCH --time=0-12:0:00
# The job command(s):


file=$1
d=$2  #number of PCs to output

/software/flashpca_x86-64 --bfile $file -d $d

