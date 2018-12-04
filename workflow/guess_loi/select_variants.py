import argparse
from subprocess import check_call

from workflow.guess_loi.checks import check_gatk


def select_variants():
    parser = argparse.ArgumentParser(description="""
        A wrapper for select variants from GATK
        """)
    parser.add_argument('output', help="output file name")
    parser.add_argument('vcf', help="vcfs of the SNPs")
    parser.add_argument('genome', help="the reference genome")
    parser.add_argument('add', nargs='*', help="additional arguments for select variants."
                                               "For example, use --select-type-to-exclude INDEL or "
                                               "--restrict-alleles-to BIALLELIC")

    args = parser.parse_args()
    gatk = check_gatk()
    run_select_variants(gatk, args.genome, args.vcf, args.additional, args.output)


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
