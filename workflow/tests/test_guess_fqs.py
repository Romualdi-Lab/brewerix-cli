from os.path import abspath
from subprocess import check_call
from tempfile import TemporaryDirectory
from typing import List

genome = abspath("./workflow/testfiles/genome.fa")

bed = abspath("./workflow/testfiles/genes_regions.bed")
snps = abspath("./workflow/testfiles/sample-contig.vcf")

fastq_file = abspath("./workflow/testfiles/single_trueSeq.fq.gz")
fastq_1 = abspath("./workflow/testfiles/paired_trueSeq_1.fastq.gz")
fastq_2 = abspath("./workflow/testfiles/paired_trueSeq_2.fastq.gz")


def guess_loi_builder(args: List, opts: List, mode: str = 'fqs'):
    cmd = ["guess-LOI", mode]
    cmd += opts
    cmd += args
    return cmd


def test_guess_single_end(hisat2_genome_index):
    with TemporaryDirectory() as wdir:
        cmd = guess_loi_builder([snps, bed, genome, hisat2_genome_index, fastq_file], [], mode='fqs')
        check_call(cmd, cwd=wdir)


def test_guess_paired_end(hisat2_genome_index):
    with TemporaryDirectory() as wdir:
        cmd = guess_loi_builder([snps, bed, genome, hisat2_genome_index, fastq_1, fastq_2], ['--paired'], mode='fqs')
        check_call(cmd, cwd=wdir)


def test_guess_bam(hisat2_bam):
    with TemporaryDirectory() as wdir:
        cmd = guess_loi_builder([snps, bed, genome, hisat2_bam], [], mode='bams')
        check_call(cmd, cwd=wdir)
