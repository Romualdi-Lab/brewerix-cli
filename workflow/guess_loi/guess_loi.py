from os.path import expanduser
from subprocess import call, DEVNULL
from sys import argv

from workflow.guess_loi.checks import check_gatk
from workflow.guess_loi.format_ase import format_ase_internal


def guess_loi():
    print(argv)
    bams = argv[3:]
    SNPs = argv[1]
    bed = argv[2]
    samples = [ bam.rstrip(".bam") for bam in bams ]

    print(bams)
    print(SNPs)
    print(bed)
    print(samples)

    gatk = check_gatk(gatk="~/local/stow/gatk-4.0.4.0/gatk")

    # check Genome Index
    # sort bam if need be
    # create RG tags, if need be

    aser_count(gatk, bams, SNPs, bed, samples)
    format_ase_cicle(samples)
    collapse_ase(samples)
    get_ASE_table()

def aser_count(gatk, bams, SNPs, bed, samples):
    # run aser count
    # format the results
    # annotate with bed file
    # create allele count info (compress the output)
    for bam, sample in zip(bams, samples):
        run_aser_on_bam(gatk, bam, SNPs, "genome-sort.fa", sample)



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
    all_keys = []
    for sample in samples:
        intermediate_file = ''.join([sample, ".aser"])

        with open(intermediate_file, 'r') as fd:
            d = {}
            for line in fd:
                key, value = line.rstrip('\n').split('\t')
                d[key] = value
            all_dictionaries.append(d)
            all_keys = set(all_keys + d.keys())
    return all_dictionaries, all_keys


def collapse(all_dictionaries, all_keys, tbl):
    for k in all_keys:
        line_values = [k]
        for dict in all_dictionaries:
            value = dict[k] if k in dict else './.'
            line_values.append(value)
        tbl.write('\t'.join(line_values) + '\n')


if __name__ == '__main__':
    guess_loi()

