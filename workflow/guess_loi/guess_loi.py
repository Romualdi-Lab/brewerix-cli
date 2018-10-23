from collections import Counter
from os.path import basename, exists
from sys import argv

from workflow.guess_loi.ase import get_ase_table, ase_table
from workflow.guess_loi.checks import check_gatk, check_file_exists
from workflow.guess_loi.common_utility import guess_sample_name, check_paired_end_nomenclature
from workflow.guess_loi.filter_count_compress_output import sort_file_by_gene_name
from workflow.guess_loi.hisat2_align import align_genome_paired_end, align_genome_single_end
from workflow.guess_loi.samtools_wrappers import samtools_filter, check_rg_tag
from workflow.guess_loi.snp_gene_association import annotate_aser_table_from_bed


def guess_loi_from_bams():
    snps = argv[1]
    bed = argv[2]
    genome = argv[3]
    bams = argv[4:]
    samples = [bam.rstrip(".bam") for bam in bams]

    gatk = check_gatk(gatk="~/local/stow/gatk-4.0.4.0/gatk")

    for bam in bams:
        check_file_exists(bam)
        if check_rg_tag(bam):
            print("error: bam with no RG tag found: " + bam)
            exit(976)

    create_ase_table_from_bams(snps, bams, bed, gatk, genome, samples)


def guess_loi_from_fqs():
    snps = argv[1]
    bed = argv[2]
    genome_idx = argv[3]
    genome = argv[4]
    mode = argv[5]
    fqs = argv[6:]

    gatk = check_gatk(gatk="~/local/stow/gatk-4.0.4.0/gatk")

    samples = []
    if mode == "paired":
        for fq in fqs:
            check_paired_end_nomenclature(fq)
            samples.append(guess_sample_name(fq, paired=True))

        samples = check_sample_paired_end(samples)

        for sample in samples:
            fq1 = sample + "_1.fq.gz"
            fq2 = sample + "_2.fq.gz"
            spl = basename(sample)
            bam_file = spl + '.bam'
            if not exists(bam_file):
                align_genome_paired_end(fq1, fq2, '4', genome_idx)
                out_file = spl + '.filter.bam'
                samtools_filter(out_file, bam_file, bed, thread_sam='4')

    elif mode == "single":
        samples = [guess_sample_name(fq) for fq in fqs]
        for sample in samples:
            fq = sample + ".fq.gz"
            spl = basename(sample)
            bam_file = spl + '.bam'
            if not exists(bam_file):
                align_genome_single_end(fq, '4', genome_idx)
                out_file = spl + '.filter.bam'
                samtools_filter(out_file, bam_file, bed, thread_sam='4')

    else:
        print("error: unrecognized mode. Choose either single or paired")
        exit(334)

    samples = [basename(s) for s in samples]
    bams = [s + '.filter.bam' for s in samples]
    create_ase_table_from_bams(snps, bams, bed, gatk, genome, samples)


def create_ase_table_from_bams(snps, bams, bed, gatk, genome, samples):
    table = ase_table(gatk, bams, snps, genome, samples)

    annotate_aser_table_from_bed(table, bed, "LOI-table.txt")
    sort_file_by_gene_name("LOI-table.txt", "LOI-table-sort.txt")
    get_ase_table("LOI-table-sort.txt")


def check_sample_paired_end(samples):
    samples_counts = Counter(samples)
    for sample in samples_counts:
        if samples_counts[sample] != 2:
            print("a sample does not have hsi pair: " + sample)
            exit(278)
    return samples_counts.keys()


# def aser_count(gatk, bams, snps, genome, samples):
#     # run aser count
#     # format the results
#     # annotate with bed file
#     # create allele count info (compress the output)
#     for bam, sample in zip(bams, samples):
#         run_aser_on_bam(gatk, bam, snps, genome, sample)


if __name__ == '__main__':
    guess_loi_from_fqs()
