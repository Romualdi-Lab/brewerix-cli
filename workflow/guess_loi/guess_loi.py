import locale
from os.path import exists, join
from re import sub
from tempfile import TemporaryDirectory
from typing import List

import pkg_resources

from workflow.guess_loi.alignments import align
from workflow.guess_loi.ase import ase_table
from workflow.guess_loi.checks import check_file_exists, check_command_availability
from workflow.guess_loi.concat_vcfs import run_concat_vcfs
from workflow.guess_loi.filter_count_compress_output import filter_useful_snps, \
    compute_overall_expression
from workflow.guess_loi.haplotype_caller_rna import run_haplotype_caller
from workflow.guess_loi.parse_args import parse_args
from workflow.guess_loi.progress import Progress
from workflow.guess_loi.samples import paired_samples, single_samples, Sample
from workflow.guess_loi.samtools import check_rg_tag, call_samtools_index
from workflow.guess_loi.select_variants import run_select_variants
from workflow.guess_loi.snp_gene_association import annotate_aser_table_from_bed, read_bed_index, create_gene2tss
from workflow.guess_loi.table import create_guess_loi_table, collapse_to_gene_info
from workflow.guess_loi.vcf_related_functions import annotate_vcf_with_heterozygous_genotype


class InputError(Exception):
    pass


def guess_loi():
    args = parse_args()
    setup_locale()
    check_command_availability(["gatk", "samtools", "bcftools", "hisat2"])

    version = pkg_resources.require("workflow-guess-loi")[0].version

    print("Workflow guess LOI version: " + version)
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
        create_ase_table_from_bams(args.snps, args.multi, args.bams, args.bed, args.genome_dict, samples, p,
                                   args.threads)


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

        # if False:
        #     lp = LineProfiler()
        #     lp_create_ase_table_from_bams = lp(create_ase_table_from_bams)
        #     lp_create_ase_table_from_bams(args.snps, args.multi, bams, args.bed, args.genome_dict, samples, p)
        #     lp.print_stats()
        # else:
        create_ase_table_from_bams(args.snps, args.multi, bams, args.bed, args.genome_dict, samples, p, args.threads)


def split_threads(threads):
    if threads == 1:
        return 1, 1
    else:
        s = threads // 2
        return threads, s


def create_annotated_lines(informative: list, overall: dict, gene_col: int, genes2tss: dict,
                           snp_id_text: str = "rs_multi", fake_pvalue: float = 1.0):
    all_genes_set = set(overall.keys())

    for line in informative:
        line_cp = line[:]
        gene = line_cp[gene_col]
        gset = {gene}
        gene_tss = genes2tss[gene]
        all_genes_set = all_genes_set.difference(gset)
        line_cp[1] = int(line_cp[1])
        line_cp.insert(gene_col+1, int(gene_tss))
        yield line_cp

    for gene in all_genes_set:
        annotation = overall[gene][0] + overall[gene][1]
        gene_tss = int(genes2tss[gene])
        collapsed_gene = collapse_to_gene_info(annotation, overall[gene][2], snp_id_text, fake_pvalue)
        collapsed_gene.insert(gene_col+1, gene_tss)
        yield collapsed_gene


def create_ase_table_from_bams(snps, multi_snps, bams, bed, genome, samples, progress, threads):
    # TODO: implement the quantification of the expression with htseq
    names = [s.name for s in samples]

    if multi_snps is None:
        vcf = snps
    else:
        vcf = resolve_multi_snps(snps, multi_snps, genome, bams, progress)

    table = ase_table(bams, vcf, genome, names, progress, threads)

    bed_idx = read_bed_index(bed)
    snp_lines = annotate_aser_table_from_bed(table, bed_idx)

    next(snp_lines)
    informative_snps = filter_useful_snps(snp_lines, gene_col=5, ratio_min=0.1)

    # The generator needs to be re-created
    snp_lines = annotate_aser_table_from_bed(table, bed_idx)
    header = next(snp_lines)
    gene_expression_estimates = compute_overall_expression(snp_lines, gene_col=5)

    gene2tss = create_gene2tss(bed)

    annotated_lines = create_annotated_lines(informative_snps, gene_expression_estimates,
                                             gene_col=5, genes2tss=gene2tss)
    header.insert(6, "TSS")

    # TODO: bypass this ordering step.
    # lines = sort_file_by_gene_name_and_position(annotated_lines)

    progress.start('Format results')
    create_guess_loi_table(annotated_lines, header, "guess_loi_table.txt")
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
