from subprocess import check_call
from sys import argv

from workflow.guess_loi.checks import check_command_availability


def hisat_build_index():
    mfasta = argv[1]  # multifasta file
    pref = argv[2]  # prefix
    build_index(mfasta, pref)


def build_index(multi_fasta, prefix="index"):
    check_command_availability(['hisat2'])
    check_call(["hisat2-build", multi_fasta, prefix])


if __name__ == '__main__':
    hisat_build_index()
