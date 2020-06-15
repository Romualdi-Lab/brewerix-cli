import locale
from multiprocessing.pool import Pool
from os import mkdir
from os.path import exists, join
from re import sub
from typing import List

import pkg_resources

from workflow.guess_loi.alignments import align
from workflow.guess_loi.ase import ase_table
from workflow.guess_loi.checks import check_file_exists, check_command_availability
from workflow.guess_loi.concat_vcfs import run_concat_vcfs
from workflow.guess_loi.filter_count_compress_output import filter_useful_snps, \
    compute_overall_expression
from workflow.guess_loi.haplotype_caller_rna import haplotype_wrapper
from workflow.guess_loi.parse_args import parse_args
from workflow.guess_loi.progress import Progress
from workflow.guess_loi.samples import paired_samples, single_samples, Sample
from workflow.guess_loi.samtools import check_rg_tag, call_samtools_index, split_bam_by_chromosomes, \
    split_vcf_by_chromosomes
from workflow.guess_loi.select_variants import run_select_variants
from workflow.guess_loi.snp_gene_association import annotate_aser_table_from_bed, read_bed_index, create_gene2tss, \
    create_gene2info
from workflow.guess_loi.table import create_guess_loi_table, collapse_to_gene_info
from workflow.guess_loi.vcf_related_functions import annotate_vcf_with_heterozygous_genotype, remove_invalid_snv_ids


class InputError(Exception):
    pass


def guess_loi():
    args = parse_args()
    setup_locale()
    check_command_availability(["gatk", "samtools", "bcftools", "hisat2"])

    version = pkg_resources.require("brewerix-cli")[0].version

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
                                   args.threads, args.gatkmem)


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

        create_ase_table_from_bams(args.snps, args.multi, bams, args.bed, args.genome_dict, samples, p,
                                   args.threads. args.gatkmem)


def split_threads(threads):
    if threads == 1:
        return 1, 1
    else:
        s = threads // 2
        return threads, s


def create_annotated_lines(informative: List, overall: dict, gene_col: int,
                           genes2tss: dict, genes2info: dict,
                           snp_id_text: str = "rs_multi", fake_pvalue: float = 1.0) -> List:

    all_genes_set = set(overall.keys())

    for line in informative:
        line_cp = line[:]
        gene = line_cp[gene_col]
        gset = {gene}

        insert_info = genes2info[gene] + [int(genes2tss[gene])]
        all_genes_set = all_genes_set.difference(gset)
        line_cp[1] = int(line_cp[1])

        line_cp = line_cp[:gene_col+1] + insert_info + line_cp[gene_col+1:]
        yield line_cp

    for gene in all_genes_set:
        annotation = overall[gene][0] + overall[gene][1]
        insert_info = genes2info[gene] + [int(genes2tss[gene])]

        collapsed_gene = collapse_to_gene_info(annotation, overall[gene][2], snp_id_text, fake_pvalue)
        collapsed_gene = collapsed_gene[:gene_col+1] + insert_info + collapsed_gene[gene_col+1:]

        yield collapsed_gene


def create_ase_table_from_bams(snps, multi_snps, bams, bed, genome, samples, progress, threads, gatkmem):
    names = [s.name for s in samples]
    bed_idx = read_bed_index(bed)

    if multi_snps is None:
        vcf = snps
    else:
        chromosomes = list(bed_idx.keys())
        vcf = resolve_multi_snps(snps, multi_snps, genome, bams, progress, chromosomes, threads)

    table = ase_table(bams, vcf, genome, names, progress, threads, gatkmem)

    snp_lines = annotate_aser_table_from_bed(table, bed_idx)

    next(snp_lines)
    informative_snps = filter_useful_snps(snp_lines, gene_col=5, ratio_min=0.1)

    # The generator needs to be re-created
    snp_lines = annotate_aser_table_from_bed(table, bed_idx)
    header = next(snp_lines)

    gene_expression_estimates = compute_overall_expression(snp_lines, gene_col=5)

    gene2tss = create_gene2tss(bed)
    gene2info = create_gene2info(bed)

    annotated_lines = create_annotated_lines(informative_snps, gene_expression_estimates,
                                             gene_col=5, genes2tss=gene2tss, genes2info=gene2info)

    header = header[:6] + ["type", "source", "TSS"] + header[6:]
    # TODO: bypass this ordering step.

    progress.start('Format results')
    create_guess_loi_table(annotated_lines, header, "brewer-table.txt")
    progress.complete()


def _resolve_multi_snps(snps: str, multi_snps: str, genome: str, bams: List[str], progress: Progress,
                        chromosomes: List[str], threads: int = 1) -> str:
    # DEPRECATED
    output = "hc-merged.vcf"
    if not exists(output):
        progress.start("Resolve multi SNPs")

        # with TemporaryDirectory() as wdir:
        if True:
            wdir = "intermediate_file"
            mkdir(wdir)

            file_list = [split_vcf_by_chromosomes(multi_snps, wdir, chromosomes)]
            for bam in bams:
                file_list.append(split_bam_by_chromosomes(bam, wdir, chromosomes))

            file_chunks = list(zip(*file_list))

        # exit(1)
        #
        # with TemporaryDirectory() as wdir:

            with Pool(threads) as pool:
                args = ((chrom_files, genome, join(wdir, str(idx) + '.vcf'))
                        for idx, chrom_files in enumerate(file_chunks))

                files = pool.map(haplotype_wrapper, args)

            called = join(wdir, "called.vcf")
            concatenated = join(wdir, "concatenated.vcf")
            called_gt = join(wdir, "called_gt.vcf")

            run_concat_vcfs([snps] + files, concatenated)
            annotate_vcf_with_heterozygous_genotype(called, called_gt, "Genotype")

            # run_haplotype_caller(multi_snps, bams, genome, called)
            # annotate_vcf_with_heterozygous_genotype(called, called_gt, "Genotype")
            # run_concat_vcfs([snps, called_gt], concatenated)
            run_select_variants(genome, concatenated, ["--select-type-to-exclude", "INDEL"], output)

        progress.complete()

    return output


def resolve_multi_snps(snps: str, multi_snps: str, genome: str, bams: List[str], progress: Progress,
                       chromosomes: List[str], threads: int = 1) -> str:

    output = "hc-merged.vcf"
    if not exists(output):
        progress.start("Resolve multi SNPs")

        # with TemporaryDirectory() as wdir:
        if True:
            wdir = "intermediate_file"
            called = join(wdir, "called.vcf")
            if not exists(wdir):
                mkdir(wdir)

            if not exists(called):
                args = ((multi_snps, bams, genome, join(wdir, chrom + '.vcf'), chrom, 1)
                        for chrom in set(chromosomes) - {"X"})

                with Pool(threads) as pool:
                    files = pool.map(haplotype_wrapper, args)
                    if "X" in chromosomes:
                        files.append(haplotype_wrapper((multi_snps, bams, genome, join(wdir, 'X.vcf'), "X", threads)))

                called = join(wdir, "called.vcf")
                run_concat_vcfs(files, called)

            concatenated = join(wdir, "concatenated.vcf")
            called_gt = join(wdir, "called_gt.vcf")

            if not exists(called_gt):
                annotate_vcf_with_heterozygous_genotype(called, called_gt, "Genotype")

            if not exists(concatenated):
                run_concat_vcfs([snps, called_gt], concatenated)
                remove_invalid_snv_ids(concatenated)

            run_select_variants(genome, concatenated, ["--select-type-to-exclude",
                                                       "INDEL", '--exclude-ids', '.'],
                                output)

        progress.complete()

    return output


if __name__ == '__main__':
    guess_loi()
