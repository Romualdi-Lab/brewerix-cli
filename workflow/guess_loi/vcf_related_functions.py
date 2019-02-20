import argparse
from typing import List


def reformat_vcf_with_het_genotype():
    parser = argparse.ArgumentParser(description="""
                Create a vcf with all the SNPs with heterozyous genotype.
                """)
    parser.add_argument('output', help="output file name")
    parser.add_argument('vcf', help="vcfs of the SNPs")
    parser.add_argument('--sampleName', help="the header name for GT column", default="gentype")

    args = parser.parse_args()

    annotate_vcf_with_heterozygous_genotype(args.vcf, args.output, args.sampleName)


def annotate_vcf_with_heterozygous_genotype(vcf, vcf_out, sampleName):
    with open(vcf, "rt") as vf, open(vcf_out, "wt") as vo:
        for line_idx, line in enumerate(vf):
            if line[:2] == "##":
                vo.write(line)
            else:
                tokens = line.rstrip().split('\t')
                if line[:1] == "#":
                    new_line = tokens[:9] + [sampleName + "\n"]  # type: List
                else:
                    new_line = tokens[:8] + ["GT", "0/1\n"]  # type: List

                vo.write('\t'.join(new_line))
