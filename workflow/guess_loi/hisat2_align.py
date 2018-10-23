from os.path import basename
from subprocess import call
from sys import argv
from tempfile import NamedTemporaryFile

from workflow.guess_loi.checks import check_file_exists
from workflow.guess_loi.common_utility import guess_sample_name
from workflow.guess_loi.samtools_wrappers import samtools_sort


def run_hisat2():
    if len(argv) == 1:
        print("run-hisat2 paired|single fq1 [fq2] threads index")
        exit(0)
    mode = argv[1]
    if mode == 'paired':
        if len(argv) != 6:
            print("Unexpected parameter")
            exit(434)

        fq1 = argv[2]
        fq2 = argv[3]
        threads = argv[4]
        index = argv[5]
        align_genome_paired_end(fq1, fq2, threads, index)

    elif mode == 'single':
        if len(argv) != 5:
            print("Unexpected parameter")
            exit(434)

        fq = argv[2]
        threads = argv[3]
        index = argv[4]
        align_genome_single_end(fq, threads, index)
    else:
        print("run-hisat2 paired|single fq1 [fq2] threads index")


def align_genome_single_end(fq, threads, index_dir_prefix, thread_sam='1'):
    check_file_exists(fq)

    sample = basename(guess_sample_name(fq))
    output_file = sample + '.bam'
    # TODO:  --rna-strandness RF evaluate strandness

    with NamedTemporaryFile(suffix=".sam") as sam:
        print(sam.name)
        cmd = ['hisat2',
               '-p', threads,
               '--rg-id', sample,
               '--rg', "LB:" + sample,
               '--rg', "PL:platform",
               '--rg', "PU:S_" + sample,
               '--rg', "SM:" + sample,
               '-x', index_dir_prefix,
               '-U', fq,
               '-S', sam.name]
        print(' '.join(cmd))
        call(cmd)

        samtools_sort(output_file, sam, thread_sam)


def align_genome_paired_end(fq1, fq2, threads, index_dir_prefix, thread_sam='1'):
    check_file_exists(fq1)
    check_file_exists(fq2)

    sample = guess_sample_name(fq1, paired=True)

    # TODO:  --rna-strandness RF evaluate strandness
    output_file = basename(sample) + '.bam'
    with NamedTemporaryFile(suffix=".sam") as sam:
        print(sam.name)
        cmd = ['hisat2',
               '-p', threads,
               '--rg-id', sample,
               '--rg', "LB:" + sample,
               '--rg', "PL:platform",
               '--rg', "PU:" + sample,
               '--rg', "SM:" + sample,
               '-x', index_dir_prefix,
               '-1', fq1,
               '-2', fq2,
               '-S', sam.name]
        print(' '.join(cmd))
        call(cmd)

        samtools_sort(output_file, sam, thread_sam)


if __name__ == '__main__':
    run_hisat2()
