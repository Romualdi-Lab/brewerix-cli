import argparse
from subprocess import check_call


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
    run_select_variants(args.genome, args.vcf, args.add, args.output)


def run_select_variants(genome, vcf, additional, output):
    # --select-type-to-exclude INDEL \

    cmd = [
              'gatk', "SelectVariants",
              "-R", genome,
              "-V", vcf,
              "-O", output,
              "--exclude-non-variants",
          ] + additional

    print(' '.join(cmd))
    check_call(cmd)


if __name__ == '__main__':
    select_variants()
