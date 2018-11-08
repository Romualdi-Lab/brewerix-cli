from collections import Counter
from os.path import basename, exists

from workflow.guess_loi.ase import create_guess_loi_table, ase_table
from workflow.guess_loi.checks import check_gatk, check_file_exists
from workflow.guess_loi.common_utility import guess_sample_name, check_paired_end_nomenclature, \
    guess_sample_name_from_bam
from workflow.guess_loi.filter_count_compress_output import sort_file_by_gene_name_and_position
from workflow.guess_loi.hisat2_align import align_genome_paired_end, align_genome_single_end
from workflow.guess_loi.parse_args import parse_args
from workflow.guess_loi.samtools_wrappers import samtools_filter, check_rg_tag
from workflow.guess_loi.snp_gene_association import annotate_aser_table_from_bed


class InputError(Exception):
    pass


def guess_loi():
    args = parse_args()
    if args.mode == 'bams':
        guess_loi_from_bams(args)
    elif args.mode == 'fqs':
        guess_loi_from_fqs(args)
    else:
        exit("unrecognized mode: %r "), args.mode


def guess_loi_from_bams(args):
    gatk = check_gatk(gatk="~/local/stow/gatk-4.0.4.0/gatk")

    for bam in args.bams:
        check_file_exists(bam)
        if check_rg_tag(bam):
            raise InputError("error: bam with no RG tag found: %r", bam)

    samples = [guess_sample_name_from_bam(bam) for bam in args.bams]
    create_ase_table_from_bams(args.snps, args.bams, args.bed, gatk, args.genome_dict, samples)


def guess_loi_from_fqs(args):
    gatk = check_gatk(gatk="~/local/stow/gatk-4.0.4.0/gatk")

    samples = []
    if args.is_paired:
        for fq in args.fqs:
            check_paired_end_nomenclature(fq)
            samples.append(guess_sample_name(fq, paired=True))

        samples = check_sample_paired_end(samples)

        for sample in samples:
            fq1 = sample + "_1.fq.gz"
            fq2 = sample + "_2.fq.gz"
            spl = basename(sample)
            bam_file = spl + '.bam'
            if not exists(bam_file):
                align_genome_paired_end(fq1, fq2, '4', args.genome_idx)
                out_file = spl + '.filter.bam'
                samtools_filter(out_file, bam_file, args.bed, thread_sam='4')

    else:
        samples = [guess_sample_name(fq) for fq in args.fqs]
        for sample in samples:
            fq = sample + ".fq.gz"
            spl = basename(sample)
            bam_file = spl + '.bam'
            if not exists(bam_file):
                align_genome_single_end(fq, '4', args.genome_idx)
                out_file = spl + '.filter.bam'
                samtools_filter(out_file, bam_file, args.bed, thread_sam='4')

    samples = [basename(s) for s in samples]
    bams = [s + '.filter.bam' for s in samples]
    create_ase_table_from_bams(args.snps, bams, args.bed, gatk, args.genome_dict, samples)


def create_ase_table_from_bams(snps, bams, bed, gatk, genome, samples):
    table = ase_table(gatk, bams, snps, genome, samples)
    annotated_lines = annotate_aser_table_from_bed(table, bed)
    head, lines = sort_file_by_gene_name_and_position(annotated_lines)
    create_guess_loi_table(lines, head, 5, "guess_loi_table.txt")


def check_sample_paired_end(samples):
    samples_counts = Counter(samples)
    for sample in samples_counts:
        if samples_counts[sample] != 2:
            raise InputError("a sample does not have its pair: %r", sample)
    return samples_counts.keys()


if __name__ == '__main__':
    guess_loi()
