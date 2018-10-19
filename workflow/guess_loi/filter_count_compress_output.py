from itertools import groupby
from operator import itemgetter
from sys import argv


def compact_snps():
    table = argv[1]
    with open(table, 'r') as fd:
        compact_snps_core(fd, gene_col=5)

def compact_snps_core(fd, gene_col=5):
    first_line = True

    for key, values in groupby(iter_snp(fd, gene_col), itemgetter(gene_col)):

        # TODO: use the iterator rather than the list?
        values = list(values)

        if first_line:
            # line = values[0]
            # print
            # '\t'.join(line)
            first_line = False
            continue

        interesting_snp = []
        overall_gene_expression = []

        for line in values:
            allelic_expressions = line[(gene_col + 1):]
            a_ratio = [ratio_allele_expression(x) for x in allelic_expressions]
            a_sum = [sum_allele_expression(x) for x in allelic_expressions]
            overall_gene_expression.append(a_sum)

            for r in a_ratio:
                if r >= 0.1:
                    interesting_snp.append(line)
                    break

        overall_gene = [sum(x) for x in zip(*overall_gene_expression)]

        if len(interesting_snp) == 0:
            line = values[0]
            refAltFake = zip(overall_gene, [0 for i in range(len(overall_gene))])
            overall_gene = [','.join([str(v) for v in x]) for x in refAltFake]
            annot = line[:(gene_col + 1)]
            annot[2] = "rs_multi"
            output = annot + overall_gene

            print('\t'.join(output)) # TODO change this to save on file
        else:
            for line in interesting_snp:
                print('\t'.join(line))

def iter_snp(fd, gene_col):
    # TODO: control if input is sorted
    for line_idx, line in enumerate(fd):
        tokens = line.rstrip('\n').split('\t')
        if len(tokens) <= gene_col:
            exit('Insufficient token number at line %d: %s' % (line_idx + 1, line))
        yield tokens


def sum_allele_expression(ae):
    if ae == "NA":
        return 0

    aes = ae.split(',')
    return sum([int(x) for x in aes])


def ratio_allele_expression(ae):
    if ae == 'NA':
        return 0

    aes = ae.split(',')
    aes = [int(x) for x in aes]

    return min(aes) / float(sum(aes))


if __name__ == '__main__':
    compact_snps()
