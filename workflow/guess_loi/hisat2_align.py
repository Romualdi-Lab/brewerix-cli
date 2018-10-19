from re import sub
from subprocess import call

from workflow.guess_loi.checks import check_file_exists


def align_genome_single_end(fq, threads, index_dir_prefix, thread_sam=1):
    check_file_exists(fq)

    sample = sub(r'\.fq(\.gz)*$', '', fq)
    output_file = sample + '.bam'
    # TODO:  --rna-strandness RF evaluate strandness

    call(['hisat2',
          '-p', threads,
          '-x', index_dir_prefix,
          '-U', fq,
          '| samtools view -u',
          '| samtools sort -l 9 -m 1G --threads', thread_sam,
          '-o', output_file])

def align_genome_paired_end(fq1, fq2, threads, index_dir_prefix, thread_sam=1):
    check_file_exists(fq1)
    check_file_exists(fq2)

    sample = sub(r'_R?1\.fq(\.gz)*$', '', fq1)

    # TODO:  --rna-strandness RF evaluate strandness
    output_file = sample + '.bam'
    call(['hisat2',
          '-p', threads,
          '-x', index_dir_prefix,
          '-1', fq1,
          '-2', fq2,
          '| samtools view -u',
          '| samtools sort -l 9 -m 1G --threads', thread_sam,
          '-o', output_file])

