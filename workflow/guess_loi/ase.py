from collections import OrderedDict
from subprocess import check_call

from workflow.guess_loi.filter_count_compress_output import compact_snps_core


def ase_table(gatk, bams, snps, genome, samples):
    ase_tables = OrderedDict([(sample, aser_count(gatk, bam, snps, genome, sample))
                              for bam, sample in zip(bams, samples)])
    return collapse_ase(ase_tables)


def aser_count(gatk, bam, vcf, genome, sample):
    step1 = sample + ".ASER.txt"
    step2 = sample + ".aser"

    check_call([
        gatk, "ASEReadCounter",
        "-I", bam,
        "-V", vcf,
        "-R", genome,
        "-O", step1,
    ])
    format_ase(step1, step2)

    return step2


def format_ase(input_, output):
    with open(input_, 'rt') as i, open(output, "wt") as o:
        format_ase_internal(i, o)


def format_ase_internal(fd, intermediate, min_coverage=5):
    for idx, line in enumerate(fd):
        tokens = line.rstrip('\n').split('\t')
        if len(tokens) < 7:
            print("less than 7 tokens at line " + str(idx + 1))

        gene_id = '_'.join(tokens[:5])
        ref = tokens[5]
        alt = tokens[6]
        if idx == 0:
            ref_alt = ','.join([ref, alt])
            intermediate.write(gene_id + '\t' + ref_alt + '\n')
        else:
            if int(ref) + int(alt) >= min_coverage:
                ref_alt = ','.join([ref, alt])
                intermediate.write(gene_id + '\t' + ref_alt + '\n')


def collapse_ase(tables):
    output = "ASER_table.txt"

    all_dictionaries, all_keys = create_keys_dictionaris(tables.values())
    with open(output, 'wt') as tbl:
        header = "chr_pos_rs_ref_alt\t" + '\t'.join(tables.keys())
        tbl.write(header + '\n')
        collapse(all_dictionaries, all_keys, tbl)

    return output


def get_ase_table(file="ASER_table.txt", gene_col=5, out_file_name="final-output-table.txt"):
    with open(file, 'r') as tbl:
        compact_snps_core(tbl, gene_col, out_file_name)


def create_keys_dictionaris(files):
    all_dictionaries = []
    all_keys = set()
    for file_ in files:
        with open(file_, 'r') as fd:
            d = {}
            for line in fd:
                key, value = line.rstrip('\n').split('\t')
                d[key] = value
            all_dictionaries.append(d)
            all_keys = all_keys | set(d.keys())
    return all_dictionaries, all_keys


def collapse(all_dictionaries, all_keys, tbl):
    for k in all_keys:
        line_values = [k]
        for dictionary in all_dictionaries:
            value = dictionary[k] if k in dictionary else 'NA'
            line_values.append(value)
        tbl.write('\t'.join(line_values) + '\n')
