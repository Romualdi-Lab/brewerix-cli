from os.path import basename
from subprocess import check_call
from sys import argv

from workflow.guess_loi.checks import check_gatk


def select_variants():
    if len(argv) == 1 or len(argv) < 3:
        print("Run like this:")
        print('%s OUTPUT VCF GENOME ADDITIONALS\n' % basename(argv[0]))
        print("For example, use --select-type-to-exclude INDEL or --restrict-alleles-to BIALLELIC")
        exit(0)

    gatk = "~/local/stow/gatk-4.0.4.0/gatk"
    output = argv[1]
    vcf = argv[2]
    genome = argv[3]
    additional = []
    if (len(argv) >3 ):
        additional = argv[4:]

    run_select_variants(gatk, genome, vcf, additional, output)


def run_select_variants(gatk, genome, vcf, additional, output):
    # --select-type-to-exclude INDEL \

    gatk = check_gatk(gatk)
    cmd = [
        gatk, "SelectVariants",
        "-R", genome,
        "-V", vcf,
        "-O", output,
        "--exclude-non-variants",
    ] + additional

    print(' '.join(cmd))
    check_call(cmd)
    

if __name__ == '__main__':
    select_variants()
