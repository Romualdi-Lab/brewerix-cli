from os.path import join, exists
from subprocess import check_call
from tempfile import TemporaryDirectory

from workflow.guess_loi.general import remove_files
from workflow.guess_loi.vcf_related_functions import safe_rename


def split_N_cigar_reads(bams, samples, genome, progress, clean=False):
    split_bams = []
    with TemporaryDirectory(dir=".") as wdir:
        for bam, sample in progress.track('SplitNCigarReads', zip(bams, samples)):
            split_bam = sample.name + ".split.bam"

            if exists(split_bam):
                print("Split bam file already found. Skip %s" % split_bam)
                split_bams.append(split_bam)
                continue

            tmp_split_bam = join(wdir, split_bam)
            check_call(["gatk", "SplitNCigarReads",
                        "-R", genome,
                        "-I", bam,
                        "-O", tmp_split_bam])

            safe_rename(tmp_split_bam, split_bam)
            split_bams.append(split_bam)

    if clean:
        remove_files(bams)

    return split_bams