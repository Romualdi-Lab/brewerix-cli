# BrewerIX-cli

## Detect loss of imprinting or X-chromosome inactivation erosion from RNAseq data

## Prerequisite

### Software (must be available on PATH)
* GATK (>= 4.0.0)
* hisat2
* samtools

### Knowledge-base
To run brewerix-cli you need to create a knowledgebase for your species of interest.
The knowledge-base must be composed as follows: 
  * a BED file with the interesting regions
  * the list of biallelic SNVs that lay on the potentially imprinted genes and sexual chromosomes
  * the list of multi-allelic SNVs that lay on the potentially imprinted genes and sexual chromosomes (optional)
  * the genome (with both genome dict and genome fa index '.fai')
  * the genome indexed for hisat2

We build a wrapper to create the knowledge-base called brewerix-prepare-knowledgebase.
Minimal requirements are 

## General Workflow

This pipeline is able to generate a table of bi-allelic SNVs starting from fastq files.

In the following lines we summarize the major steps:
* collect fastqs
* generate the bams filtered by the region of interest
* Haplotype caller for multi-allelic SNVs (optional)
  * merge bi-allelic SNVs with the bi-allelic batch (optional)
* process with ASEReadCounted (GATK) on bi-allelic SNVs
* merge the sample in the global matrix
* annotation and additional filtering

The final output is a table that can be analyzed by brewerIX app.

## Quick start

The simplest situation is to perform the analysis of single end samples\
with the standard pipeline, i.e. using a pre-built set of bi-allelic SNVs. 

The call should look like this one

guess-LOI fqs snps-biallelic.vcf \ <br />
&ensp; genes-and-sex-regions.bed \ <br />
&ensp; genome.fa hisat2-index/genome *.fastq.gz

