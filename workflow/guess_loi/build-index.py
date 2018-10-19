from subprocess import call

from workflow.guess_loi.checks import check_hisat2_installation


def build_index(multi_fasta, prefix="index"):
    check_hisat2_installation()
    call(["hisat2-build", multi_fasta, prefix])
