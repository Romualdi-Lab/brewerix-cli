from collections import Counter
from subprocess import call
from sys import argv

from workflow.guess_loi.checks import check_gatk
from workflow.guess_loi.common_utility import guess_sample_name, check_paired_end_nomenclature
from workflow.guess_loi.filter_count_compress_output import compact_snps_core
from workflow.guess_loi.format_ase import format_ase_internal
from workflow.guess_loi.hisat2_align import align_genome_paired_end, align_genome_single_end
from workflow.guess_loi.samtools_wrappers import samtools_filter, check_rg_tag


def guess_loi_from_bams():
    SNPs = argv[1]
    bed = argv[2]
    genome = argv[3]
    bams = argv[4:]
    samples = [ bam.rstrip(".bam") for bam in bams ]

    print(bams)
    print(SNPs)
    print(bed)
    print(samples)

    gatk = check_gatk(gatk="~/local/stow/gatk-4.0.4.0/gatk")

    for bam in bams:
        if check_rg_tag(bam):
            print("error: bam with no RG tag found: " + bam)
            exit(976)

    aser_count(gatk, bams, SNPs, genome, samples)
    format_ase_cicle(samples)
    collapse_ase(samples)
    get_ASE_table()

def guess_loi_from_fqs():
    SNPs = argv[1]
    bed = argv[2]
    genome = argv[3]
    mode = argv[4]
    fqs = argv[5:]

    gatk = check_gatk(gatk="~/local/stow/gatk-4.0.4.0/gatk")

    if mode == "paired":
        samples = []
        for fq in fqs:
            check_paired_end_nomenclature(fq)
            samples.append(guess_sample_name(fq, paired=True))

        samples = check_sample_paired_end(samples)

        for sample in samples:
            fq1 = sample + "_1.fq.gz"
            fq2 = sample + "_2.fq.gz"
            align_genome_paired_end(fq1, fq2, 4, genome)
            out_file = sample + '.filter.bam'
            input = sample + '.bam'
            samtools_filter(out_file, input, bed, thread_sam=4)

    elif mode == "single":
        samples = [ guess_sample_name(fq) for fq in fqs ]
        for sample in samples:
            fq = sample + ".fq.gz"
            align_genome_single_end(fq, 4, genome)
            out_file = sample + '.filter.bam'
            input = sample + '.bam'
            samtools_filter(out_file, input, bed, thread_sam=4)

    else:
        print("error: unrecognized mode. Choose either single or paired")
        exit(334)

    bams = [s + '.filter.bam' for s in samples]
    aser_count(gatk, bams, SNPs, genome, samples)
    format_ase_cicle(samples)
    collapse_ase(samples)
    get_ASE_table()




def check_sample_paired_end(samples):
    samples_counts = Counter(samples)
    for sample in samples_counts:
        if samples_counts[sample] != 2:
            print("a sample does not have hsi pair: " + sample)
            exit(278)
    return samples_counts.keys()


def aser_count(gatk, bams, SNPs, genome, samples):
    # run aser count
    # format the results
    # annotate with bed file
    # create allele count info (compress the output)
    for bam, sample in zip(bams, samples):
        run_aser_on_bam(gatk, bam, SNPs, genome, sample)



def run_aser_on_bam(gatk, bam, vcf, genome, sample):
    call([gatk, "ASEReadCounter",
          "-I", bam,
          "-V", vcf,
          "-R", genome,
          "-O", ''.join([ sample,".ASER.txt"])])


def get_ASE_table(file="ASER_table.txt"):
    with open(file, 'r') as tbl:
        compact_snps_core(tbl, gene_col=5)


def format_ase_cicle(samples):
    for sample in samples:
        intermediate_file = ''.join([sample, ".ASER.txt"])
        with open(intermediate_file, 'r') as fd:
            intermediate_out = ''.join([sample, ".aser"])
            with open(intermediate_out, "w") as intermediate:
                format_ase_internal(fd, intermediate)


def collapse_ase(samples):
    all_dictionaries, all_keys = create_keys_dictionaris(samples)
    intermediate_table_file = "ASER_table.txt"
    with open(intermediate_table_file, 'w') as tbl:
        collapse(all_dictionaries, all_keys, tbl)


def create_keys_dictionaris(samples):
    all_dictionaries = []
    all_keys = set()
    for sample in samples:
        intermediate_file = ''.join([sample, ".aser"])

        with open(intermediate_file, 'r') as fd:
            d = {}
            for line in fd:
                key, value = line.rstrip('\n').split('\t')
                d[key] = value
            all_dictionaries.append(d)
            all_keys = set(all_keys + set(d.keys()))
    return all_dictionaries, all_keys


def collapse(all_dictionaries, all_keys, tbl):
    for k in all_keys:
        line_values = [k]
        for dict in all_dictionaries:
            value = dict[k] if k in dict else './.'
            line_values.append(value)
        tbl.write('\t'.join(line_values) + '\n')


if __name__ == '__main__':
    guess_loi_from_bams()

