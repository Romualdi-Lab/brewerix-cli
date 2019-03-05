import argparse
from os.path import join
from subprocess import check_call
from tempfile import TemporaryDirectory
from typing import List

from workflow.guess_loi.checks import check_command_availability


def concat_vcfs():
    parser = argparse.ArgumentParser(description="""
            Concatenate VCFs
            """)
    parser.add_argument('output', help="output file name")
    parser.add_argument('vcfs', nargs='+', help="vcfs of the SNPs")

    args = parser.parse_args()

    run_concat_vcfs(args.vcfs, args.output)


def run_concat_vcfs(files: List, output: str):
    check_command_availability(['bcftools', 'bgzip', 'tabix'])

    with TemporaryDirectory() as wdir:
        gzipped_files = []

        for file in files:
            gfile = join(wdir, file + ".gz")
            gzipped_files.append(gfile)

            with open(gfile, "wb") as dest:
                check_call(["bgzip", "-c", file], stdout=dest)
            check_call(["tabix", "-p", "vcf", gfile])

        cmd = [
                  'bcftools', "concat",
                  "-a", "-D",
                  "-O", "v",
                  "-o", output,
              ] + gzipped_files
        print(' '.join(cmd))
        check_call(cmd)
