from os.path import exists

from workflow.guess_loi.alignments import align
from workflow.guess_loi.ase import create_guess_loi_table, ase_table
from workflow.guess_loi.checks import check_gatk, check_file_exists
from workflow.guess_loi.filter_count_compress_output import sort_file_by_gene_name_and_position
from workflow.guess_loi.parse_args import parse_args
from workflow.guess_loi.progress import Progress
from workflow.guess_loi.samples import paired_samples, single_samples
from workflow.guess_loi.samtools import check_rg_tag, call_samtools_index
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
    raise RuntimeError('broken')
    gatk = check_gatk(gatk=args.gatk)

    for bam in args.bams:
        check_file_exists(bam)
        if check_rg_tag(bam):
            raise InputError("error: bam with no RG tag found: %r", bam)

    samples = [guess_sample_name_from_bam(bam) for bam in args.bams]
    with Progress(args.progress) as p:
        create_ase_table_from_bams(args.snps, args.bams, args.bed, gatk, args.genome_dict, samples, p)


def guess_loi_from_fqs(args):
    gatk = check_gatk(gatk=args.gatk)
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
            call_samtools_index(bam)

        create_ase_table_from_bams(args.snps, bams, args.bed, gatk, args.genome_dict, samples, p)


def split_threads(threads):
    # TODO: better heuristic needed
    if threads == 1:
        return 1, 1
    else:
        s = threads // 2
        return threads, s


def create_ase_table_from_bams(snps, bams, bed, gatk, genome, samples, progress):
    # TODO: implement the quantification of the expression with htseq
    names = [s.name for s in samples]
    table = ase_table(gatk, bams, snps, genome, names, progress)
    annotated_lines = annotate_aser_table_from_bed(table, bed)
    head, lines = sort_file_by_gene_name_and_position(annotated_lines)

    progress.start('Format results')
    create_guess_loi_table(lines, head, 5, "guess_loi_table.txt")
    progress.complete()


if __name__ == '__main__':
    guess_loi()
