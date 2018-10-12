# GUESS LOI

## Infer from RNAseq data the loss of imprinting of imprinted genes.

## Available species
* Human hg38
* Mouse m38

## Prerequisite

* Create a BED file with the interesting regions
* Create the list of SNPs that lay on the potentially imprinted genes and sexual cromosome
    - Human / Mouse
    - use gatk utility IndexFeatureFile
    - Use Paolo's utility extract biallelic.sh [PICARD?]
* Mind the pseudoautosomal regions
    - Human / Mouse

## General Pipiline Workflow

* collect the bams
* check genome index presence and version
* sort bams if not sorted
* create the RG tags if not present
* process with ASER (RNAseq-preprocessing-and-ASER-gatk4.sh)
    - inputs: potentially imprinted SNPs and BAMS
* format_imprinting_ASERcount.sh + extract_info-from-id.sh to create the results.
* annotate-snos with BED file of imprinted genes
* filter-snp-allele-counts.py
