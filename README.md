# RGStraP
[![DOI](https://zenodo.org/badge/443880532.svg)](https://zenodo.org/badge/latestdoi/443880532) [![Snakemake](https://img.shields.io/badge/snakemake-≥6.15.5-brightgreen.svg?style=flat)](https://snakemake.readthedocs.io)

**RGStraP** (***R***NA-seq-based ***G***enetic ***Stra***tification ***P***Cs) is a bioinformatics pipeline for calculating Principal Components (PCs) showing genetic stratification from RNA-seq data. The pipeline mainly utilizes the variant calling capabilities of [GATK4](https://gatk.broadinstitute.org/hc/en-us/articles/360035531192-RNAseq-short-variant-discovery-SNPs-Indels-) and the principal component analysis (PCA) of [FlashPCA2](https://github.com/gabraham/flashpca). The pipeline was built using [snakemake](https://snakemake.github.io/).

![RNAvc_Figure1](https://user-images.githubusercontent.com/30294080/174513895-7c18c769-2488-4e40-be00-0ce9011a9708.png)


## Contact
Muhamad Fachrul, [mfachrul@student.unimelb.edu.au](mailto:mfachrul@student.unimelb.edu.au?subject=[GitHub]%20RGStraP)

[![Twitter Follow](https://img.shields.io/twitter/follow/f_azr?style=social&logo=twitter)](https://twitter.com/f_azr)

## Citation
Coming soon

## Requirements
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
- Modify the `config/config.yaml` file according to where the necessary files are in your system. Variables to modify include:
  - Path to a file containing list of ONLY the first pair of paired-end fastq samples to be analyzed.
  - Path to metadata file (required for adding read-group information with GATK).
    - Has to be a tab-delimited file with 6 columns and no header, with the first column containing BAM file locations with the format `2_mapped/[FILENAME]_Aligned.sortedByCoord.out.bam` and the next five columns representing read-group ID, platform, sample name, library, and platform unit, respectively.
    - More info [here](https://gatk.broadinstitute.org/hc/en-us/articles/360035890671-Read-groups), example here.
  - Path to directory of [reference genome index generated by STAR](https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf).
  - Path to reference genome fasta file.
  - Path to indel files (for GATK's BaseRecalibrator).
  - Path to **flashpca**.
- Please adjust the 'dupedist' value according to your sequencing platform in the `scripts/clumpify_OpDup.sh` file (recommendations included within the script).
- Test the pipeline by performing a dry-run.
```
snakemake -n
```
- Running the pipeline on a cluster using a workload manager / job scheduler, such as [slurm](https://slurm.schedmd.com/documentation.html), is highly recommended. An example of a snakemake profile to run it on slurm is included.
  - Please modify the partition name in `slurm/config.yaml` file accordingly.
  - You can also modify the maximum number of jobs to be run at once in the `slurm/config.yaml` file.
```
# To run pipeline on slurm
snakemake --profile slurm
```
### Running the pipeline from VCF file (lite version)
RGStraP can also be used to capture RG-PCs from existing VCF files via the lite version.
- Make sure to modify the `config/lite_config.yaml` file accordingly.
```
# To run lite pipeline on slurm
snakemake -s lite_Snakefile --cores 2
```
## License
Apache 2.0 License
