## RGStraP

import itertools
import pandas as pd
import numpy as np

configfile: "config/config.yaml"

# Get sample name and initial directory    
fastq_list = pd.read_csv(config["sample_fq_list"], header=None, index_col=False, squeeze = True)
SAMPLES = fastq_list.str.rsplit('/', 1, expand=True).rename(columns={0: "dir",1: "name"})
SAMPLES["ext"] = SAMPLES["name"].str.split('.', 1, expand=True)[1]
SAMPLES["name"] = SAMPLES["name"].str.rsplit('1.f', 1, expand=True)[0]

# create folder for cluster use
onstart:
    print("##### Creating profile pipeline #####\n") 
    print("\t Creating jobs output folder...\n")
    shell("mkdir -p jobs/")

# Paths to 

rule all:
    input:
        fqc_zip = expand("fastqc/{name}1_fastqc.zip", name = SAMPLES["name"]),
        fqc_html = expand("fastqc/{name}1_fastqc.html", name = SAMPLES["name"]),
	trimmed_1 = expand("trimmed/{name}1_val_1.fq.gz", name = SAMPLES["name"]),
        trimmed_2 = expand("trimmed/{name}2_val_2.fq.gz", name = SAMPLES["name"]),
        clumped_1 = expand("no_OpDup/{name}1_val_1_clumped.fq.gz", name = SAMPLES["name"]),
        clumped_2 = expand("no_OpDup/{name}2_val_2_clumped.fq.gz", name = SAMPLES["name"]),
        pass1_bam = expand("1_mapped/{name}Aligned.sortedByCoord.out.bam", name = SAMPLES["name"]),
        pass1_sj = expand("1_mapped/{name}SJ.out.tab", name = SAMPLES["name"]),
        pass1_finLog = expand("1_mapped/{name}Log.final.out", name = SAMPLES["name"]),
        pass2_bam = expand("2_mapped/{name}Aligned.sortedByCoord.out.bam", name = SAMPLES["name"]),
        pass2_sj = expand("2_mapped/{name}SJ.out.tab", name = SAMPLES["name"]),
        pass2_finLog = expand("2_mapped/{name}Log.final.out", name = SAMPLES["name"]),
        withRG = expand("addRG/{name}withRG.bam", name = SAMPLES["name"]),
        filt = expand("addRG/{name}withRG_filtered.bam", name = SAMPLES["name"]),
        nodup = expand("nodup/{name}withRG_filtered_marked.bam", name = SAMPLES["name"]),
        ncigar = expand("NCIGAR/{name}withRG_filtered_marked_splitN.bam", name = SAMPLES["name"]),
        recaltab = expand("bqsr/{name}recal.table", name = SAMPLES["name"]),
        finalbam = expand("bqsr/{name}final.bam", name = SAMPLES["name"]),
        firvcf = expand("vars/{name}first.g.vcf.gz", name = SAMPLES["name"]),
        sj_list = "pass1_sj.list",
        gvcf_list = "Merged.gvcf.list",
        genotyped = "VCF/Genotyped_raw.vcf.gz",
        varFil = "VCF/Genotyped_filterKeep.vcf.gz",
        varSel = "VCF/Genotyped_filterOut.vcf.gz",
        bim_pre = "VCF/Genotyped_filterOut.bim",
        bim_clean = "VCF/Genotyped_filterOut_clean.bim",
        bim_maf = "VCF/Genotyped_filterOut_clean_maf005.bim",
        hm3_list = "VCF/overlap_maf005_hapmap3_RNA_SNPs.list",
        hm3_bim = "VCF/Genotyped_filterOut_clean_maf005_hm3.bim",
        ld_bim = "VCF/Genotyped_filterOut_clean_maf005_hm3_ld005.bim",
        pcs = "VCF/RG-StraP_pcs.txt",
        eigenvalues = "VCF/RG-StraP_eigenvalues.txt"
    output: touch("all.done")

# Run QC before mapping (fastqc and trimming adapters)
rule pre_QC:
    input:
        p1_raw = SAMPLES["dir"][0]+"/{name}1."+SAMPLES["ext"][0]
    output:
        fqc_zip = "fastqc/{name}1_fastqc.zip",
        fqc_html = "fastqc/{name}1_fastqc.html",
	trimmed_1 = "trimmed/{name}1_val_1.fq.gz",
        trimmed_2 = "trimmed/{name}2_val_2.fq.gz"

    resources:
        mem_mb=6000
        
    threads: 4
        
    shell:
        """
        sh scripts/fastqc.sh {input.p1_raw}
	sh scripts/trim.sh {input.p1_raw}
        """

# Run BBMap's Clumpify to remove optical duplicates pre-mapping
## Please adjust the 'dupedist' value according to your sequencing platform in the scripts/clumpify_OpDup.sh file
## (recommendations included within the script).
rule clumpify:
    input:
        trimmed_1 = rules.pre_QC.output.trimmed_1
    
    output:
        clumped_1 = "no_OpDup/{name}1_val_1_clumped.fq.gz",
        clumped_2 = "no_OpDup/{name}2_val_2_clumped.fq.gz"

    resources:
        mem_mb=60000  
 
    threads: 4
       
    shell:
        """
	sh scripts/clumpify_OpDup.sh {input.trimmed_1}
	"""

# Run first-pass mapping of RNA-seq files to reference genome with STAR
rule onePass_STAR:
    input:
        clumped_1 = rules.clumpify.output.clumped_1
        
    output:
        pass1_bam = "1_mapped/{name}Aligned.sortedByCoord.out.bam",
        pass1_sj = "1_mapped/{name}SJ.out.tab",
        pass1_finLog = "1_mapped/{name}Log.final.out"

    resources:
        mem_mb=100000      
 
    params:
        genome = config["STAR_genome"]

    threads: 8
    
    shell:
        """
        sh scripts/1pass_STAR.sh {input.clumped_1} {params.genome}
        """

# Get all splice junction files from first-pass mapping before second-pass
rule getSJs:
    input:
        expand("1_mapped/{name}SJ.out.tab", name = SAMPLES["name"])

    output:
        sj_list = "pass1_sj.list"

    shell:
        """
        echo "{input}" > {output.sj_list}
        """

# Run second-pass mapping of RNA-seq files to reference genome with STAR
rule twoPass_STAR:
    input:
        clumped_1 = rules.clumpify.output.clumped_1,
        sj_list = rules.getSJs.output.sj_list

    output:
        pass2_bam = "2_mapped/{name}Aligned.sortedByCoord.out.bam",
        pass2_sj = "2_mapped/{name}SJ.out.tab",
        pass2_finLog = "2_mapped/{name}Log.final.out"

    resources:
        mem_mb=100000

    params:
        genome = config["STAR_genome"]

    threads: 8

    shell:
        """
        sh scripts/2pass_STAR.sh {input.clumped_1} {params.genome} {input.sj_list}
        """

# Add read group information to BAM files, which are necessary during GATK's variant calling process.
rule addRG:
    input:
        pass2_bam = rules.twoPass_STAR.output.pass2_bam

    output:
        withRG = "addRG/{name}withRG.bam"

    resources:
        mem_mb=20000

    params:
        meta = config["meta"]
 
    threads: 8
    
    shell:
        """
        sh scripts/addRG.sh {input.pass2_bam} {params.meta} scripts/PICARD_addRG.sh
        """

# Filter BAM files based on quality with samtools
rule filterBAM:
    input:
        withRG = rules.addRG.output.withRG

    output:
        filt = "addRG/{name}withRG_filtered.bam"

    resources:
        mem_mb=5000

    threads: 8

    shell:
        """
        samtools view -bq 20 {input.withRG} > {output.filt}
        """

# Sort BAM files, then mark and remove duplicates
rule sortMarkDup:
    input:
        filt = rules.filterBAM.output.filt

    output:
        nodup = "nodup/{name}withRG_filtered_marked.bam"

    resources:
        mem_mb=30000

    threads: 8

    shell:
        """
        sh scripts/sortMarkDup.sh {input.filt}
        """

# Splits reads that contain Ns in their cigar string
rule splitN:
    input:
        nodup = rules.sortMarkDup.output.nodup

    output:
        ncigar = "NCIGAR/{name}withRG_filtered_marked_splitN.bam"

    resources:
        mem_mb=20000

    params:
        genome_fa = config["genome_fa"]

    threads: 8

    shell:
        """
        sh scripts/splitncigar.sh {input.nodup} {params.genome_fa}
        """

# Perform base recalibration (get the recalibration table)
rule baseRecalib:
    input:
        nodup = rules.splitN.output.ncigar

    output:
        recaltab = "bqsr/{name}recal.table"

    resources:
        mem_mb=20000

    params:
        genome_fa = config["genome_fa"],
        indel1 = config["indel1"],
        indel2 = config["indel2"]

    threads: 8

    shell:
        """
        sh scripts/base_recalibrator.sh {input.nodup} {params.genome_fa} {params.indel1} {params.indel2}
        """

# Apply base recalibration
rule applybqsr:
    input:
        nodup = rules.splitN.output.ncigar,
        recaltab = rules.baseRecalib.output.recaltab

    output:
        finalbam = "bqsr/{name}final.bam"

    resources:
        mem_mb=40000

    params:
        genome_fa = config["genome_fa"]

    threads: 8

    shell:
        """
        sh scripts/apply_bqsr.sh {input.nodup} {params.genome_fa} {input.recaltab}
        """

# Call germline SNPs and indels via local re-assembly of haplotypes
rule haploCall:
    input:
        finalbam = rules.applybqsr.output.finalbam

    output:
        "vars/{name}first.g.vcf.gz"

    resources:
        mem_mb=50000

    params:
        genome_fa = config["genome_fa"]

    threads: 8

    shell:
        """
        sh scripts/haploCall.sh {input.finalbam} {params.genome_fa}
        """

# Get list of GVCF files for merging
rule getListGeno:
    input:
        expand("vars/{name}first.g.vcf.gz", name = SAMPLES["name"])

    output:
        gvcf_list = "Merged.gvcf.list"

    shell:
        """
        echo "{input}" | sed 's/vars/-V vars/g' > {output.gvcf_list}
        """

# Merge GVCFs and call genotype
rule genotype:
    input:
        gvcf_list = rules.getListGeno.output.gvcf_list

    output:
        genotyped = "VCF/Genotyped_raw.vcf.gz"

    resources:
        mem_mb=80000,
        runtime="48h"

    params:
        genome_fa = config["genome_fa"],
        interval = config["interval"]

    threads: 8

    shell:
        """
        sh scripts/genotypegvcfs.sh {input.gvcf_list} {params.genome_fa} {params.interval}
        """

# Initial filtering of variants called
rule filGen:
    input:
        genotyped = rules.genotype.output.genotyped

    output:
        varFil = "VCF/Genotyped_filterKeep.vcf.gz"

    resources:
        mem_mb=20000

    params:
        genome_fa = config["genome_fa"]

    threads: 8

    shell:
        """
        sh scripts/variantFil.sh {input.genotyped} {params.genome_fa}
        """

# Select variants based on previous filtering
rule selGen:
    input:
        varFil = rules.filGen.output.varFil

    output:
        varSel = "VCF/Genotyped_filterOut.vcf.gz"

    resources:
        mem_mb=20000

    params:
        genome_fa = config["genome_fa"]

    threads: 8

    shell:
        """
        sh scripts/variantSel.sh {input.varFil} {params.genome_fa}
        """

# Clean up VCF file
rule plink_clean:
    input:
        varSel = rules.selGen.output.varSel

    output:
        bim_pre = "VCF/Genotyped_filterOut.bim",
        bim_clean = "VCF/Genotyped_filterOut_clean.bim"

    threads: 1

    shell:
        """
        plink --vcf {input.varSel} --make-bed --out VCF/Genotyped_filterOut
        sh scripts/idtochrpos.sh {output.bim_pre}
        Rscript scripts/clean_data.R VCF/Genotyped_filterOut VCF/Genotyped_filterOut_clean
        """

# Filter based on MAF<0.05
rule qc_plink:
    input:
        bim_clean = rules.plink_clean.output.bim_clean
        
    output:
        bim_maf = "VCF/Genotyped_filterOut_clean_maf005.bim"

    threads: 1

    shell:
        """
        sh scripts/qc_all_basic.sh VCF/Genotyped_filterOut_clean VCF/Genotyped_filterOut_clean_maf005 0.05 0 0.1 0.2
        """

# Compare with HapMap3 reads and take overlaps
rule hapmap3:
    input:
        bim_maf = rules.qc_plink.output.bim_maf

    output:
        hm3_list = "VCF/overlap_maf005_hapmap3_RNA_SNPs.list",
        hm3_bim = "VCF/Genotyped_filterOut_clean_maf005_hm3.bim"

    params:
        hapmap3 = config["hapmap3"]

    threads: 1

    shell:
        """
        cut -f1,4,5,6 {input.bim_maf} > temp_ChrID.txt
        echo -e "CHROM\tPOS\tREF\tALT" | cat - temp_ChrID.txt > RNA_ChrID.txt
        rm ./temp_ChrID.txt

        Rscript scripts/hapmap3_overlap.R RNA_ChrID.txt {params.hapmap3}

        plink --bfile ./VCF/Genotyped_filterOut_clean_maf005 --extract {output.hm3_list} --make-bed --out ./VCF/Genotyped_filterOut_clean_maf005_hm3
        """

# Filter based on ld<0.05
rule ld:
    input:
        hm3_bim = rules.hapmap3.output.hm3_bim

    output:
        ld_bim = "VCF/Genotyped_filterOut_clean_maf005_hm3_ld005.bim"

    params:
        ld_reg = config["ld_regions"]

    threads: 1

    shell:
        """
        sh scripts/ld_thinning_pruning.sh VCF/Genotyped_filterOut_clean_maf005_hm3 VCF/Genotyped_filterOut_clean_maf005_hm3_ld005 1000 50 0.05 {params.ld_reg}
        """

# Generate genetic PCs
rule pca:
    input:
        ld_bim = rules.ld.output.ld_bim

    output:
        pcs = "VCF/RG-StraP_pcs.txt",
        eigenvalues = "VCF/RG-StraP_eigenvalues.txt"

    params:
        flashpca = config["flashpca"]

    shell:
        """
        {params.flashpca} --bfile VCF/Genotyped_filterOut_clean_maf005_hm3_ld005 -d 10 --outpc VCF/RG-StraP_pcs.txt --outpve VCF/RG-StraP_pve.txt --outvec VCF/RG-StraP_eigenvectors.txt --outval VCF/RG-StraP_eigenvalues.txt
        """
