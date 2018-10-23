from itertools import groupby
from operator import itemgetter
from sys import argv


def compact_snps():
    table = argv[1]
    gene_col = int(argv[2])
    with open(table, 'r') as fd:
        compact_snps_core(fd, gene_col, "/dev/stdout")

def compact_snps_core(fd, gene_col=5, out_file_name="final-output-table.txt"):
    with open(out_file_name, "w") as out:
        first_line = True

        for key, values in groupby(iter_snp(fd, gene_col), itemgetter(gene_col)):

            # TODO: use the iterator rather than the list?
            values = list(values)

            if first_line:
                out.write('\t'.join(values[0]) + '\n')
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

                # print('\t'.join(output)) # TODO change this to save on file
                out.write('\t'.join(output) + '\n')
            else:
                for line in interesting_snp:
                    # print('\t'.join(line))
                    out.write('\t'.join(line) + '\n')


def sort_file_by_gene_name(in_file, out_file):
    with open(in_file, "r") as fd, open(out_file, "w") as out:
        d = {}
        for idx, line in enumerate(fd):
            if idx == 0:
                out.write(line)
                continue

            tks = line.split('\t')
            key = tks[5]

            if key in d:
                d[key].append(line)
            else:
                d[key] = [line]

        for k in sorted(d.keys()):
            str = ''.join(d[k])
            out.write(str)


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
