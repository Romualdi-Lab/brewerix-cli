import argparse
from logging import warning
from os import rename
from os.path import join, basename, dirname
from tempfile import TemporaryDirectory
from typing import List

from pysam.libcbcf import VariantFile


def reformat_vcf_with_het_genotype():
    parser = argparse.ArgumentParser(description="""
                Create a vcf with all the SNPs with heterozyous genotype.
                """)
    parser.add_argument('output', help="output file name")
    parser.add_argument('vcf', help="vcfs of the SNPs")
    parser.add_argument('--sampleName', help="the header name for GT column", default="gentype")

    args = parser.parse_args()

    annotate_vcf_with_heterozygous_genotype(args.vcf, args.output, args.sampleName)


def annotate_vcf_with_heterozygous_genotype(vcf, vcf_out, sample_name):
    with open(vcf, "rt") as vf, open(vcf_out, "wt") as vo:
        for line_idx, line in enumerate(vf):
            if line[:2] == "##":
                vo.write(line)
            else:
                tokens = line.rstrip().split('\t')
                if line[:1] == "#":
                    new_line = tokens[:9] + [sample_name + "\n"]  # type: List
                else:
                    new_line = tokens[:8] + ["GT", "0/1\n"]  # type: List

                vo.write('\t'.join(new_line))


def exclude_invalid_snv_ids():
    parser = argparse.ArgumentParser(description="""
                Exclude ids containing "_" from vcf. 
                """)
    parser.add_argument('vcf_file', help="vcf file name")
    args = parser.parse_args()

    remove_invalid_snv_ids(args.vcf_file)


def remove_invalid_snv_ids(vcf_file):
    vcf = VariantFile(vcf_file)
    vcf_filename = basename(vcf_file)
    with TemporaryDirectory(dir=".") as wdir:
        temp_file = join(wdir, vcf_filename)
        with VariantFile(temp_file, 'w', header=vcf.header) as out:
            for line in vcf:
                snv_id = line.id
                if "_" in snv_id:
                    warning("Invalid SNV id found will be filtered out: %s" % snv_id)
                    continue
                out.write(line)

        safe_rename(temp_file, vcf_file)


def safe_rename(source, target):
    with TemporaryDirectory(dir=dirname(target)) as temp:
        intermediate = join(temp, 'obj')
        rename(source, intermediate)
        rename(intermediate, target)
