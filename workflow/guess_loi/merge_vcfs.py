from os.path import basename, join
from subprocess import check_call
from sys import argv
from tempfile import TemporaryDirectory


def merge_vcfs():
    if len(argv) == 1 or len(argv) < 2:
        print("Run like this:")
        print('%s OUTPUT VCFs\n' % basename(argv[0]))
        exit(0)

    output = argv[1]
    files = argv[2:]

    run_merge_vcfs(files, output)


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
