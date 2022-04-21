# RGStraP
**RGStraP** (***R***NA-seq-based ***G***enetic ***Stra***tification ***P***Cs) is a bioinformatics pipeline for calculating Principal Components (PCs) showing genetic stratification from RNA-seq data. The pipeline mainly utilizes the variant calling capabilities of [GATK4](https://gatk.broadinstitute.org/hc/en-us/articles/360035531192-RNAseq-short-variant-discovery-SNPs-Indels-) and the principal component analysis (PCA) of [FlashPCA2](https://github.com/gabraham/flashpca). The pipeline was built using [snakemake](https://snakemake.github.io/).

<img src='https://user-images.githubusercontent.com/30294080/156269248-866ae75a-5ac2-4643-a443-c56a8286ecd9.png' width='700'>

## Contact
Muhamad Fachrul, [mfachrul@student.unimelb.edu.au](mailto:mfachrul@student.unimelb.edu.au?subject=[GitHub]%20RGStraP)

[![Twitter Follow](https://img.shields.io/twitter/follow/f_azr?style=social&logo=twitter)](https://twitter.com/f_azr)

## Requirements
[![Snakemake](https://img.shields.io/badge/snakemake-≥6.15.5-brightgreen.svg?style=flat)](https://snakemake.readthedocs.io)

Most of the dependencies (including [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) v0.11.8, [Trim galore](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/) v0.6.0, [BBMap](https://github.com/BioInfoTools/BBMap) (for Clumpify.sh), [STAR](https://github.com/alexdobin/STAR) v2.7.10a, [Picard](https://broadinstitute.github.io/picard/) v2.24.0, [Samtools](http://www.htslib.org/) v1.8, [GATK4](https://gatk.broadinstitute.org/hc/en-us/articles/360036194592-Getting-started-with-GATK4) v4.0.6.0, and [PLINK 1.9](https://www.cog-genomics.org/plink/) v1.90b6.16) are included in the setup.

Please install [FlashPCA](https://github.com/gabraham/flashpca) v2.0 from source.

## How to use
### Installling Conda and snakemake
- Install a Conda-based Python3 distribution such as [Miniconda](https://docs.conda.io/en/latest/miniconda.html) or [Mambaforge](https://github.com/conda-forge/miniforge#mambaforge). In this case, we will use the latter as an example.
```
curl -L https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-x86_64.sh -o Mambaforge-Linux-x86_64.sh
bash Mambaforge-Linux-x86_64.sh
```
- Create and/or navigate to the directory in which you want the analysis of your project to take place, then clone this repository.
```
git clone https://github.com/fachrulm/RGStraP
```
- Change into the **RGStraP** directory, and create a Conda environment to run the pipeline.
```
cd RGStraP

# Activate Conda environment
conda activate base

# Create RGStraP environment
mamba env create --name RGStraP --file environment.yaml
```
- Activate the **RGStraP** environment. This environment needs to be active everytime you want to use the pipeline.
```
conda activate RGStraP

# To deactivate the environment
conda deactivate
```
### Running the pipeline
- Fill config file with paths to your file

- Option for slurm
