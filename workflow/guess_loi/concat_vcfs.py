import argparse
from os.path import join
from subprocess import check_call
from tempfile import TemporaryDirectory

from workflow.guess_loi.checks import check_bcftools


def concat_vcfs():
    parser = argparse.ArgumentParser(description="""
            Concatenate VCFs
            """)
    parser.add_argument('output', help="output file name")
    parser.add_argument('vcfs', nargs='+', help="vcfs of the SNPs")

    args = parser.parse_args()

    run_concat_vcfs(args.vcfs, args.output)


def run_concat_vcfs(files, output):
    bcftools = check_bcftools("bcftools")

    with TemporaryDirectory() as wdir:
        gzipped_files = []

        for file in files:
            gfile = join(wdir, file + ".gz")
            gzipped_files.append(gfile)

            with open(gfile, "wb") as dest:
                check_call(["bgzip", "-c", file], stdout=dest)
                check_call(["tabix", "-p", "vcf", gfile])

        cmd = [
                  bcftools, "merge",
                  "-a", "-D",
                  "-O", "v",
                  "-o", output,
              ] + gzipped_files
        print(' '.join(cmd))
        check_call(cmd)
