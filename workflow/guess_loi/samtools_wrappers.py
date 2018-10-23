from subprocess import call, PIPE, check_output
from tempfile import NamedTemporaryFile

from re import search


def samtools_sort(output_file, sam, thread_sam):
    with NamedTemporaryFile(suffix=".bam") as bam:
        call(["samtools", "view", "-u", "-o", bam.name, sam.name])

        sam_cmd = ["samtools", "sort", "-l", "9", "-m", "1G", "-@", thread_sam, "-o", output_file, bam.name]
        call(sam_cmd)


def samtools_sort(output_file, sam, thread_sam):
    with NamedTemporaryFile(suffix=".bam") as bam:
        call(["samtools", "view", "-u", "-o", bam.name, sam.name])

        sam_cmd = ["samtools", "sort", "-l", "9", "-m", "1G", "-@", thread_sam, "-o", output_file, bam.name]
        call(sam_cmd)


def samtools_filter(output_file, input, bed, thread_sam=4):
    call(["samtools", "view", "-@", thread_sam, "-b", "-L", bed, "-o", output_file, input])
    samtools_index(output_file)


def samtools_index(input):
    call(["samtools", "index", input])


def check_rg_tag(input):
    out = check_output(["samtools", "view", "-H", input], universal_newlines=True)
    for line in out.split("\n"):
        if not search(r'^@RG', line):
            break
        return True

    return False
