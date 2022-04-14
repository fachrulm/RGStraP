# RGStraP
**RGStraP** (***R***NA-seq-based ***G***enetic ***Stra***tification ***P***Cs) is a bioinformatics pipeline for calculating Principal Components (PCs) showing genetic stratification from RNA-seq data. The pipeline mainly utilizes the variant calling capabilities of [GATK4](https://gatk.broadinstitute.org/hc/en-us/articles/360035531192-RNAseq-short-variant-discovery-SNPs-Indels-) and the principal component analysis (PCA) of [FlashPCA2](https://github.com/gabraham/flashpca).

<img src='https://user-images.githubusercontent.com/30294080/156269248-866ae75a-5ac2-4643-a443-c56a8286ecd9.png' width='700'>

## Contact
Muhamad Fachrul, [mfachrul@student.unimelb.edu.au](mailto:mfachrul@student.unimelb.edu.au?subject=[GitHub]%20RGStraP)

[![Twitter Follow](https://img.shields.io/twitter/follow/f_azr?style=social&logo=twitter)](https://twitter.com/f_azr)

## Directory Contents
config: Configuration files necessary for the snakemake pipeline to run (the main one being config.yaml, which needs to be edited according to your input files).

scripts: Script files for different parts of the pipeline.

## Requirements
[![Snakemake](https://img.shields.io/badge/snakemake-≥5.6.0-brightgreen.svg?style=flat)](https://snakemake.readthedocs.io)

- [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) v0.11.8
- [Trim galore](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/) v0.6.0
- [BBMap](https://github.com/BioInfoTools/BBMap) (for Clumpify.sh) 
- [STAR](https://github.com/alexdobin/STAR) v2.7.10a
- [Picard](https://broadinstitute.github.io/picard/) v2.24.0
- [Samtools](http://www.htslib.org/) v1.8
- [GATK4](https://gatk.broadinstitute.org/hc/en-us/articles/360036194592-Getting-started-with-GATK4) v4.0.6.0
- [PLINK 1.9](https://www.cog-genomics.org/plink/) v1.90b6.16
- [FlashPCA](https://github.com/gabraham/flashpca) v2.0

