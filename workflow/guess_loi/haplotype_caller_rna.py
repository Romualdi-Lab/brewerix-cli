from argparse import ArgumentParser
from subprocess import check_call
from typing import List


def haplotype_caller():
    parser = ArgumentParser(description="""
                A wrapper for GATK HaplotypeCaller

                Be aware that your bams should be produced with RNAseq aligner.
                """)

    parser.add_argument('output', help="output file name")
    parser.add_argument('vcf', help="vcfs of the SNPs")
    parser.add_argument('genome', help="the reference genome")
    parser.add_argument('bams', nargs='+', help="bam files")

    args = parser.parse_args()
    run_haplotype_caller(args.vcf, args.bams, args.genome, args.output)


def haplotype_wrapper(args):
    vcf, bams, genome, haplotype_file, chrom, threads = args
    run_haplotype_caller(vcf, bams, genome, haplotype_file, chrom, threads)
    return haplotype_file


def run_haplotype_caller(alleles_vcf: str, bams: List[str], genome: str, haplotype_file: str, chrom: str, threads: int = 1):

    # TODO implement discovery mode

    bams_input = []
    for bam in bams:
        bams_input += ["-I", bam]

    cmd = [
        'gatk', "HaplotypeCaller",
        "--genotyping-mode", "GENOTYPE_GIVEN_ALLELES",
        "--max-alternate-alleles", "1",
        "-stand-call-conf", "1",
        "--native-pair-hmm-threads", str(threads),
        "-L", chrom + ":1+",
        "--alleles", alleles_vcf,
        "--dbsnp", alleles_vcf,
        ] + bams_input + [
        "-R", genome,
        "-O", haplotype_file,
    ]

    print(' '.join(cmd))
    check_call(cmd)


if __name__ == '__main__':
    haplotype_caller()
