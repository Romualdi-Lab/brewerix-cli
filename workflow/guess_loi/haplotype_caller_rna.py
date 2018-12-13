from os.path import basename
from subprocess import check_call
from sys import argv
from typing import List


def haplotype_caller():
    if len(argv) == 1 or len(argv) < 5:
        print("Run like this:")
        print('%s OUTPUT ALLELES_VCF GENOME BAMS...\n' % basename(argv[0]))
        exit(0)

    haplotype_file = argv[1]
    alleles_vcf = argv[2]
    genome = argv[3]
    bams = argv[4:]

    run_haplotype_caller(alleles_vcf, bams, genome, haplotype_file)


def run_haplotype_caller(alleles_vcf: str, bams: List[str], genome: str, haplotype_file: str):

    # TODO implement discovery mode

    bams_input = []
    for bam in bams:
        bams_input += ["-I", bam]

    cmd = [
        'gatk', "HaplotypeCaller",
        "--genotyping-mode", "GENOTYPE_GIVEN_ALLELES",
        "--max-alternate-alleles", "1",
        "-stand-call-conf", "1",
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
