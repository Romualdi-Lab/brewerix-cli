import re
from builtins import str, int

from argparse import ArgumentParser
from subprocess import check_call, check_output, STDOUT
from typing import List

from workflow.guess_loi.checks import vcf_index_exits


def haplotype_caller():
    parser = ArgumentParser(description="""
                A wrapper for GATK HaplotypeCaller

                Be aware that your bams should be produced with RNAseq aligner.
                """)

    parser.add_argument('output', help="output file name")
    parser.add_argument('vcf', help="vcfs of the SNPs")
    parser.add_argument('genome', help="the reference genome")
    parser.add_argument('bams', nargs='+', help="bam files")
    parser.add_argument('chrom', help="a chromosome to process")

    args = parser.parse_args()
    run_haplotype_caller(args.vcf, args.bams, args.genome, args.output, args.chrom)


def haplotype_wrapper(args):
    vcf, bams, genome, haplotype_file, chrom, threads = args
    if not vcf_index_exits(haplotype_file + '.idx'):
        run_haplotype_caller(vcf, bams, genome, haplotype_file, chrom, threads)
    return haplotype_file


def run_haplotype_caller(alleles_vcf: str, bams: List[str], genome: str, haplotype_file: str,
                         chrom: str, threads: int = 1):

    # TODO implement discovery mode

    bams_input = []
    for bam in bams:
        bams_input += ["-I", bam]

    "The Genome Analysis Toolkit (GATK) v4.1.4.0"

    version_lines = check_output(['gatk', 'HaplotypeCaller', '--version'], stderr=STDOUT, universal_newlines=True)

    raw_version = version_lines.split('\n')[0]
    m = re.match(r'Version:(\d\.\d\.\d+\.\d+)', raw_version)
    if m is None:
        m = re.match(r'.+\(GATK\) v(\d\.\d\.\d+\.\d+)', raw_version)
    version = m.group(1)

    cmd = ['gatk', "HaplotypeCaller"]
    if version <= '4.0.11.0':
        cmd += ["--genotyping-mode", "GENOTYPE_GIVEN_ALLELES"]

    cmd += ["--max-alternate-alleles", "1",
            "-stand-call-conf", "1",
            "--native-pair-hmm-threads", str(threads),
            "-L", chrom + ":1+",
            "--alleles", alleles_vcf,
            "--dbsnp", alleles_vcf,
            ] + bams_input + ["-R", genome, "-O", haplotype_file]

    print(' '.join(cmd))
    check_call(cmd)


if __name__ == '__main__':
    haplotype_caller()
