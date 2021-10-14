from workflow.guess_loi.general import remove_files
from workflow.guess_loi.vcf_related_functions import safe_rename
from subprocess import check_call
from os.path import join, exists
from os import environ
from tempfile import TemporaryDirectory


def mark_duplicates(bams, samples, progress, clean=False):
    bams_dedup = []
    with TemporaryDirectory(dir=".") as wdir:
        for bam, sample in progress.track('Mark duplicates', zip(bams, samples)):
            out_matrix = sample.name + "_dupStats_matrix.txt"
            out_bam = sample.name + "_dedup.bam"
            if exists(out_bam):
                print("Mark duplicates already performed. Skip %s" % out_bam)
                bams_dedup.append(out_bam)
                continue

            tmp_bam = join(wdir, sample.name + ".bam")
            tmp_matrix = join(wdir, out_matrix)
            my_env = dict(environ)
            my_env.update("USE_JDK_DEFLATER": 'true', 'USE_JDK_INFLATER': 'true')

            check_call(["picard", "MarkDuplicates",
                        "I=", bam,
                        "O=", tmp_bam,
                        "M=", tmp_matrix],
                       env=my_env)
            safe_rename(tmp_bam, out_bam)
            safe_rename(tmp_matrix, out_matrix)

            bams_dedup.append(out_bam)

    if clean:
        remove_files(bams)

    return bams_dedup