from os.path import basename
from sys import argv
from typing import List


def haplotype_caller_dna():
    if len(argv) == 1 or len(argv) < 5:
        print("Run like this:")
        print('%s OUTPUT ANNOTATED_SNP GENOME BAMS...\n' % basename(argv[0]))
        print('Be aware that your bams should be produced with DNAseq aligner such as BWA.\n')
        exit(0)

    output = argv[1]
    db_snp = argv[2]
    genome = argv[3]
    bams = argv[4:]

    run_haplotype_caller_dna(db_snp, bams, genome, output)


def run_haplotype_caller_dna(db_snp: str, bams: List[str], genome: str, output: str):

    bams_input = []
    for bam in bams:
        bams_input += ["-I", bam]

    # Evaluate if using -L targets.interval_list increase speed

    cmd = [
        'gatk', "HaplotypeCaller",
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
