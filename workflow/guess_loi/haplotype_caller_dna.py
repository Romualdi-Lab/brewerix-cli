from argparse import ArgumentParser
from typing import List


def haplotype_caller_dna():
    parser = ArgumentParser(description="""
            A wrapper for GATK HaplotypeCaller
            
            Be aware that your bams should be produced with DNAseq aligner such as BWA.
            """)
    parser.add_argument('output', help="output file name")
    parser.add_argument('vcf', help="vcfs of the SNPs")
    parser.add_argument('genome', help="the reference genome")
    parser.add_argument('bams', nargs='+', help="bam files")

    args = parser.parse_args()

    run_haplotype_caller_dna(args.vcf, args.bams, args.genome, args.output)


def run_haplotype_caller_dna(db_snp: str, bams: List[str], genome: str, output: str):

    bams_input = []
    for bam in bams:
        bams_input += ["-I", bam]

    # Evaluate if using -L targets.interval_list increase speed

    cmd = [
        'gatk', "HaplotypeCaller",
        '--genotyping-mode', 'DISCOVERY',
        '--max-alternate-alleles', '1',
        "-R", genome,
    ] + bams_input + [
        "--dbsnp", db_snp,
        "-stand_call_conf",  "30",
        "-O", output,
    ]

    print(' '.join(cmd))
    # check_call(cmd)


if __name__ == '__main__':
    haplotype_caller_dna()
