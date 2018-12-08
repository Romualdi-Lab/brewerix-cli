from re import search
from subprocess import PIPE, check_output, check_call, Popen


def call_samtools_index(bam):
    check_call(["samtools", "index", bam])


def check_rg_tag(bam_input):
    out = check_output(["samtools", "view", "-H", bam_input], universal_newlines=True)
    for line in out.split("\n"):
        if not search(r'^@RG', line):
            break
        return True

    return False


def samtools_view(stdin, bed):
    args = [
        'samtools',
        "view",
        "-u",
        "-L",
        bed
    ]

    return Popen(args, stdin=stdin, stdout=PIPE)


def samtools_sort(stdin, output, samtools_threads):
    args = [
        'samtools',
        'sort',
        '-l', '9',
        '-m', '1G',
        '-@', str(samtools_threads),
        '-o', output
    ]

    return Popen(args, stdin=stdin)
