# RGStraP
**RGStraP** (***R***NA-seq-based ***G***enetic ***Stra***tification ***P***Cs) is a bioinformatics pipeline for calculating Principal Components (PCs) showing genetic stratification from RNA-seq data. The pipeline mainly utilizes the variant calling capabilities of [GATK4](https://gatk.broadinstitute.org/hc/en-us/articles/360035531192-RNAseq-short-variant-discovery-SNPs-Indels-) and the principal component analysis (PCA) of [FlashPCA](https://github.com/gabraham/flashpca).

![RNAvc_Figure1](https://user-images.githubusercontent.com/30294080/156269248-866ae75a-5ac2-4643-a443-c56a8286ecd9.png)

## Requirements

- [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- [Trim galore](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/)
- [BBMap](https://github.com/BioInfoTools/BBMap) (for Clumpify.sh) 
- [STAR](https://github.com/alexdobin/STAR)
- [Picard](https://broadinstitute.github.io/picard/)
- [Samtools](http://www.htslib.org/)
- [GATK4](https://gatk.broadinstitute.org/hc/en-us/articles/360036194592-Getting-started-with-GATK4)
- [PLINK 1.9](https://www.cog-genomics.org/plink/)
- [FlashPCA](https://github.com/gabraham/flashpca)
