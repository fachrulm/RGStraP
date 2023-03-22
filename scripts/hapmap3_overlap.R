#!/usr/bin/Rscript

#mfachrul date: 20/01/2022

########################################################################################
# Load and install R package dependencies
########################################################################################

# Load / install dependencies
load_or_install <- function(pkg) {
  if(!require(pkg, character.only = TRUE, quietly=TRUE)) {
    install.packages(pkg, repos="https://cloud.r-project.org", verbose=FALSE, quiet=TRUE)
    library(pkg, character.only = TRUE)
  }
}
load_or_install("data.table")

args = commandArgs(trailingOnly=TRUE)

if(length(args) != 2) {
  stop("Use: Rscript hapmap3_overlap.R bim_file hapmap3_chrPos", call.=FALSE)
}

run_command <- function(command) {
  system(command)
}

## get allele pair to account for strand differences
getComp <- function(ALe) {
  paste(ifelse(ALe == "A", "T",
               ifelse(ALe == "T", "A",
                      ifelse(ALe == "G", "C", "G"))))
}

### Get ChrPos overlaps between MAF-filtered RNA-seq-based variants and HapMap3 list
RNA_ChrID <- read.table(args[1], header=TRUE, sep="\t")
RNA_ChrID$ChrID <- paste(as.character(RNA_ChrID$CHROM),RNA_ChrID$POS, sep = ":")
RNA_ChrID$PosAL <- paste(as.character(RNA_ChrID$CHROM),RNA_ChrID$POS,
			 RNA_ChrID$REF, RNA_ChrID$ALT, sep = ":")
RNA_ChrID$swPosAL <- paste(as.character(RNA_ChrID$CHROM),RNA_ChrID$POS,
                           RNA_ChrID$ALT, RNA_ChrID$REF, sep = ":")
RNA_ChrID$flPosAL <- paste(as.character(RNA_ChrID$CHROM),RNA_ChrID$POS,
			   getComp(RNA_ChrID$REF), getComp(RNA_ChrID$ALT), sep = ":")

hapmap3_SNPs <- read.table(args[2], header=TRUE, sep="\t")
hapmap3_SNPs$ChrID <- paste(as.character(hapmap3_SNPs$CHROM),hapmap3_SNPs$POS, sep = ":")
hapmap3_SNPs$PosAL <- paste(as.character(hapmap3_SNPs$CHROM),hapmap3_SNPs$POS,
                            hapmap3_SNPs$REF, hapmap3_SNPs$ALT, sep = ":")

RNA_hapmap_ChrID <- RNA_ChrID$ChrID[(RNA_ChrID$PosAL %in% hapmap3_SNPs$PosAL | RNA_ChrID$swPosAL %in% hapmap3_SNPs$PosAL | RNA_ChrID$flPosAL %in% hapmap3_SNPs$PosAL)]

write.table(RNA_hapmap_ChrID, file = "./VCF/overlap_maf005_hapmap3_RNA_SNPs.list",
	    quote = FALSE, row.names=FALSE, col.names=FALSE)
