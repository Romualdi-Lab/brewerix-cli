import locale
from os.path import exists, join
from re import sub
from tempfile import TemporaryDirectory
from typing import List

from workflow.guess_loi.alignments import align
from workflow.guess_loi.ase import create_guess_loi_table, ase_table
from workflow.guess_loi.checks import check_file_exists, check_command_availability
from workflow.guess_loi.concat_vcfs import run_concat_vcfs
from workflow.guess_loi.filter_count_compress_output import sort_file_by_gene_name_and_position
from workflow.guess_loi.haplotype_caller_rna import run_haplotype_caller
from workflow.guess_loi.parse_args import parse_args
from workflow.guess_loi.progress import Progress
from workflow.guess_loi.samples import paired_samples, single_samples, Sample
from workflow.guess_loi.samtools import check_rg_tag, call_samtools_index
from workflow.guess_loi.select_variants import run_select_variants
from workflow.guess_loi.snp_gene_association import annotate_aser_table_from_bed
from workflow.guess_loi.vcf_related_functions import annotate_vcf_with_heterozygous_genotype


class InputError(Exception):
    pass


def guess_loi():
    args = parse_args()
    setup_locale()
    check_command_availability(["gatk", "samtools", "bcftools", "hisat2"])

    if args.mode == 'bams':
        guess_loi_from_bams(args)
    elif args.mode == 'fqs':
        guess_loi_from_fqs(args)
    else:
        exit("unrecognized mode: %r" % args.mode)


def setup_locale():
    locale.setlocale(locale.LC_NUMERIC, 'POSIX')


def guess_loi_from_bams(args):
    # raise RuntimeError('broken')

    samples = []
    for bam in args.bams:
        check_file_exists(bam)
        if check_rg_tag(bam):
            raise InputError("error: bam with no RG tag found: %r", bam)
        samples.append(Sample(sub(r'\.bam$', '', bam), [], bam))

    with Progress(args.progress) as p:
        create_ase_table_from_bams(args.snps, args.multi, args.bams, args.bed, args.genome_dict, samples, p)


def guess_loi_from_fqs(args):
    hisat_threads, samtools_threads = split_threads(args.threads)

    samples = (paired_samples if args.is_paired else single_samples)(args.fqs)
    bams = []

    with Progress(args.progress) as p:
        for sample in p.track('Alignment', samples):
            bam = '%s.bam' % sample.name
            if not exists(bam):
                align(sample, args.genome_idx, args.bed, bam, hisat_threads, samtools_threads)
            bams.append(bam)

        for bam in p.track('Index generation', bams):
            if not exists(bam + '.bai'):
                call_samtools_index(bam)

        create_ase_table_from_bams(args.snps, args.multi, bams, args.bed, args.genome_dict, samples, p)


def split_threads(threads):
    if threads == 1:
        return 1, 1
    else:
        s = threads // 2
        return threads, s


def create_ase_table_from_bams(snps, multi_snps, bams, bed, genome, samples, progress):
    # TODO: implement the quantification of the expression with htseq
    names = [s.name for s in samples]

    if multi_snps is None:
        vcf = snps
    else:
        vcf = resolve_multi_snps(snps, multi_snps, genome, bams, progress)

    table = ase_table(bams, vcf, genome, names, progress)
    annotated_lines = annotate_aser_table_from_bed(table, bed)
    head, lines = sort_file_by_gene_name_and_position(annotated_lines)

    progress.start('Format results')
    create_guess_loi_table(lines, head, 5, "guess_loi_table.txt")
    progress.complete()


def resolve_multi_snps(snps: str, multi_snps: str, genome: str, bams: List[str], progress: Progress) -> str:
    output = "hc-merged.vcf"

    if not exists("hc-merged.vcf"):
        progress.start("Resolve multi SNPs")
        with TemporaryDirectory() as wdir:
            called = join(wdir, "called.vcf")
            concatenated = join(wdir, "concatenated.vcf")
            called_gt = join(wdir, "called_gt.vcf")
            run_haplotype_caller(multi_snps, bams, genome, called)
            annotate_vcf_with_heterozygous_genotype(called, called_gt, "gentype")
            run_concat_vcfs([snps, called_gt], concatenated)
            run_select_variants(genome, concatenated, ["--select-type-to-exclude", "INDEL"], output)

        progress.complete()

    return output


if __name__ == '__main__':
    guess_loi()
