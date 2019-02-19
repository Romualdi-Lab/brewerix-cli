from argparse import ArgumentParser
from os.path import join
from subprocess import check_call
from tempfile import TemporaryDirectory


def merge_vcfs():
    parser = ArgumentParser(description="""
                Merge VCFs into a single VCF
                """)
    parser.add_argument('output', help="output file name")
    parser.add_argument('vcfs', nargs='+', help="vcf files")

    args = parser.parse_args()

    run_merge_vcfs(args.vcfs, args.output)


def run_merge_vcfs(files, output):

    with TemporaryDirectory() as wdir:
        gzipped_files = []

        for file in files:
            gfile = join(wdir, file + ".gz")
            gzipped_files.append(gfile)

            with open(gfile, "wb") as dest:
                check_call(["bgzip", "-c", file], stdout=dest)
            check_call(["tabix", "-p", "vcf", gfile])

        cmd = [
                  'bcftools', "merge",
                  "-m", "none",
                  "-O", "v",
                  "-o", output,
              ] + gzipped_files
        print(' '.join(cmd))
        check_call(cmd)
