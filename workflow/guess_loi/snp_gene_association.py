import re
from itertools import groupby
from operator import itemgetter
from sys import argv


def snp2gene():
    """
    takes 2 argument:
        file (list of ids in 1st column or VCF) with SNPs) argv1
        bed_file (interesting gene in bed format) argv[2]
    """
    vcf_file = argv[1]
    bed_file = argv[2]
    snp2gene_association_from_bed(vcf_file, bed_file)


def snp2gene_association_from_vcf(vcf_file):
    with open(vcf_file, 'r') as fd:
        for idx, line in enumerate(fd):
            if line[0] == '#':
                continue

            tokens = line.rstrip('\n').split('\t')
            if tokens[7].find('GENEINFO=') == -1:
                print("No gene information fuond at line " + str(idx+1))
                exit(559)
            base_id = '_'.join(tokens[:6])
            p = re.compile(r'GENEINFO=([^;]+);')
            gene_infos = p.findall(tokens[7])

            for gene_info in gene_infos[0].split('|'):
                entrez, symbol = gene_info.split(':')
                print(base_id + '\t' + symbol + '\t' + entrez)


def snp2gene_association_from_bed(vcf_file, bed_file):
    bedidx = create_bed_index(bed_file)

    with open(vcf_file, 'r') as fd:
        for idx, line in enumerate(fd):
            if line[0] == '#':
                continue
            chr, pos, id, ref, alt, *infos= line.rstrip('\n').split('\t')
            snp_internal_id = '_'.join([chr, pos, id, ref, alt])

            if chr in bedidx:
                for start, stop, gene, _strand in bedidx[chr]:
                    if start <= pos and stop >= pos:

                        print(snp_internal_id + '\t' + gene)


def annotate():
    aserTable = argv[1]
    bed_file = argv[2]
    output_file = argv[3]
    annotate_aserTable_from_bed(aserTable, bed_file, output_file)


def annotate_aserTable_from_bed(aserTable, bed_file, output_table):
    # the first column is the id. The rest are the values.
    bedidx = create_bed_index(bed_file)

    with open(aserTable, 'r') as fd, open(output_table, 'w') as out:
        for idx, line in enumerate(fd):
            tks = line.rstrip('\n').split('\t')
            chr, pos, id, ref, alt = tks[0].split("_")
            snp_internal_id_expand = '\t'.join([chr, pos, id, ref, alt])
            values = '\t'.join(tks[1:])

            if (idx == 0):
                out.write(snp_internal_id_expand + '\t' + "symbol" + '\t' + values + '\n')
                continue

            if chr in bedidx:
                for start, stop, gene, _strand in bedidx[chr]:
                    if start <= pos and stop >= pos:
                        out.write(snp_internal_id_expand + '\t' + gene + '\t' + values + '\n')


def create_bed_index(bed_file):
    with open(bed_file, 'r') as fd:
        d = {}
        for id, grp in groupby(read_line(fd), itemgetter(0)):
            grp = [ elem[1] for elem in grp ]
            d[str(id)] = grp
        return d


def read_line(fd):
    pre_id = None
    length_t = None

    for lineidx, line in enumerate(fd):
        tokens = line.rstrip().split('\t')
        if length_t is not None and length_t != len(tokens):
            exit("Malformed input: incorrect number of columns at line %s" % (lineidx + 1))
        length_t = len(tokens)

        if pre_id is not None and tokens[0] < pre_id:
            exit("Malformed input: lexicographically sorted on col 1 at line %s" % (lineidx + 1))
        pre_id = tokens[0]

        key = tokens[0]
        info = tokens[1:]

        yield key, info

if __name__ == '__main__':
    annotate_aserTable_from_bed()
