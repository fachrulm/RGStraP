#!/usr/bin/Rscript

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
  stop("Use: Rscript 1-clean_data.R input_file output_file", call.=FALSE)
}

run_command <- function(command) {
  system(command)
}

# Function for flipping the strand of an allele.
# Uses a series of gsub calls to replace A's with T's,
# G's with C's, and vice-versa. Also works for alleles
# with more than one nucleotide (e.g. indels).
flip_strand <- function(x) {
  # Swap each letter for a dummy, we need this intermediate
  # step so we can distinguish between alleles when swapping.
  # E.g if we did A -> T then T -> A we'd end up with all A's
  # and no T's. instead we do A -> V -> T and T -> X -> A.
	x <- gsub("A", "V", x)
	x <- gsub("T", "X", x)
	x <- gsub("C", "Y", x)
	x <- gsub("G", "Z", x)
	x <- gsub("V", "T", x)
	x <- gsub("X", "A", x)
	x <- gsub("Y", "G", x)
	x <- gsub("Z", "C", x)
	return(x)
}

{
	# Load bim file
	print("Loading bim file")
	stats_fp <- paste0(args[1], ".bim")
	stats <- fread(stats_fp)
	# Name columns appropriately and extract
	tryCatch({
		setnames(stats, names(stats)[1], "chr")
		setnames(stats, names(stats)[2], "rsid")
		setnames(stats, names(stats)[4], "pos")
		setnames(stats, names(stats)[5], "A1")
		setnames(stats, names(stats)[6], "A2")
	})
	stats <- stats[, .(chr, rsid, pos, A1, A2)]
	stats$rsid <- paste0(stats$chr,":",stats$pos)

	print("Computing duplicates")
	# Drop duplicates and multiallelic variants - check both
	# by rsID and chromosome and position. I've found some cases in
	# the past where the same rsID can occur at different positions.
	dups_by_pos <- stats[,.N,by=.(rsid)][N > 1]
	stats <- stats[!dups_by_pos, on = .(rsid)]
  
	print("Removing non-SNPs")
	# Remove variants which aren't SNPS
	stats <- stats[nchar(A1) == 1 & nchar(A2) == 1]

	# Remove SNPS that are not ACGT
	bases <- c('A','C','T','G')
	stats <- stats[A1 %in% bases & A2 %in% bases]

	# Remove variants that have palindromic/strand ambiguous alleles.
	stats <- stats[A1 != flip_strand(A2)]

	# Write out list of rsIDs for plink to keep
	fwrite(stats[,.(rsid)], quote=FALSE, col.names=FALSE, file="qced_vars.txt")
}

print("Running plink")
#Extract SNPs  (maybe add maf)
cmd <- "plink1.9"
cmd <- paste(cmd, sprintf("--bfile %s", args[1]))
cmd <- paste(cmd, "--chr 1-22 --allow-extra-chr --extract qced_vars.txt")
cmd <- paste(cmd, sprintf("--make-bed --out %s", args[2]))
# Run command
system(cmd, wait=TRUE)

system("rm qced_vars.txt", wait=TRUE)

