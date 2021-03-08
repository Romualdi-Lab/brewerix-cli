from os.path import join, exists
from subprocess import check_call
from tempfile import TemporaryDirectory

from workflow.guess_loi.general import remove_files
from workflow.guess_loi.vcf_related_functions import safe_rename


def recalibration(bams, samples, genome, snvs, progress, clean=False):
    recal_bams = []

    with TemporaryDirectory(dir=".") as wdir:
        for bam, sample in progress.track('BaseRecalibrator', zip(bams, samples)):

            out_bam = sample.name + ".bqsr.bam"
            out_table = sample.name + "_recal_data.table"
            out_plot = sample.name + "_analyze_covariates.pdf"

            if exists(out_bam) and exists(out_table) and (out_plot):
                print("Base recalibration file already found. Skip %s" % out_bam)
                recal_bams.append(out_bam)
                continue

            tmp_bam = join(wdir, sample.name + ".bam")
            tmp_table = join(wdir, sample.name + "_recal_data.table")
            tmp_plot = join(wdir, sample.name + "_analyze_covariates.pdf")

            check_call(["gatk", "BaseRecalibrator",
                        "-I", bam,
                        "-R", genome,
                        "--known-sites", snvs,
                        "-O", tmp_table])

            check_call(["gatk", "ApplyBQSR",
                        "-I", bam,
                        "-R", genome,
                        "--bqsr-recal-file", tmp_table,
                        "-O", tmp_bam])

            check_call(["gatk", "AnalyzeCovariates",
                        "-bqsr", tmp_table,
                        "-plots", tmp_plot])

            safe_rename(tmp_bam, out_bam)
            safe_rename(tmp_table, out_table)
            safe_rename(tmp_plot, out_plot)

            recal_bams.append(out_bam)

    if clean:
        remove_files(bams)

    return recal_bams


def base_recalibration(bams, samples, genome, snvs, progress):
    recal_bams = []
    with TemporaryDirectory(dir=".") as wdir:
        for bam, sample in progress.track('BaseRecalibrator', zip(bams, samples)):
            split_bam = sample.name + ".split.bam"
            out_bam = sample.name + ".bqsr.bam"
            out_table = sample.name + "_recal_data.table"
            out_plot = sample.name + "_analyze_covariates.pdf"

            if exists(split_bam):
                print("Base recalibration file already found. Skip %s") % split_bam
                continue

            tmp_split_bam = join(wdir, sample.name + ".split.bam")
            check_call(["gatk", "SplitNCigarReads",
                        "-R", genome,
                        "-I", bam,
                        "-O", tmp_split_bam])

            if exists(out_bam) and exists(out_table) and (out_plot):
                print("Base recalibration file already found. Skip %s") % out_bam
                recal_bams.append(out_bam)
                continue

            tmp_bam = join(wdir, sample.name + ".bam")
            tmp_table = join(wdir, sample.name + "_recal_data.table")
            tmp_plot = join(wdir, sample.name + "_analyze_covariates.pdf")

            check_call(["gatk", "BaseRecalibrator",
                        "-I", tmp_split_bam,
                        "-R", genome,
                        "--known-sites", snvs,
                        "-O", tmp_table])
            check_call(["gatk", "ApplyBQSR",
                        "-I", tmp_split_bam,
                        "-R", genome,
                        "--bqsr-recal-file", tmp_table,
                        "-O", tmp_bam])
            check_call(["gatk", "AnalyzeCovariates",
                        "-bqsr", tmp_table,
                        "-plots", tmp_plot])

            safe_rename(tmp_split_bam, split_bam)
            safe_rename(tmp_bam, out_bam)
            safe_rename(tmp_table, out_table)
            safe_rename(tmp_plot, out_plot)

            recal_bams.append(out_bam)

    return recal_bams
