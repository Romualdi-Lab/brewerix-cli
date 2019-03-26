from os import wait, kill
from signal import SIGTERM
from subprocess import Popen, CalledProcessError, PIPE

from workflow.guess_loi.checks import check_file_exists
from workflow.guess_loi.samtools import samtools_view, samtools_sort


def align(sample, genome_index, bed, output, hisat_threads, samtools_threads):
    for fastq in sample.fastqs:
        check_file_exists(fastq)

    hisat = hisat2(sample.name, genome_index, sample.fastqs, hisat_threads)
    view = samtools_view(hisat.stdout, bed)
    sort = samtools_sort(view.stdout, output, samtools_threads)
    wait_all('hisat2', [hisat, view, sort])


def hisat2(name, genome_index, fastqs, hisat_threads):
    hisat_args = [
        'hisat2',
        '-p', str(hisat_threads),
        '--rg-id', name,
        '--rg', "LB:" + name,
        '--rg', "PL:platform",
        '--rg', "PU:" + name,
        '--rg', "SM:" + name,
        '-x', genome_index,
    ]

    assert 0 < len(fastqs) < 3
    if len(fastqs) == 1:
        hisat_args += ['-U', fastqs[0]]
    else:
        hisat_args += [
            '-1', fastqs[0],
            '-2', fastqs[1],
        ]

    return Popen(hisat_args, stdout=PIPE)


def wait_all(name, procs):
    pids = {p.pid for p in procs}

    while len(pids):
        pid, status = wait()

        assert pid in pids
        pids.remove(pid)

        if status != 0:
            kill_all(pids)

            # TODO: report the correct exit code
            raise CalledProcessError(1, name)


def kill_all(pids):
    for pid in pids:
        kill(pid, SIGTERM)

    while len(pids):
        pid = wait()[0]

        assert pid in pids
        pids.remove(pid)
