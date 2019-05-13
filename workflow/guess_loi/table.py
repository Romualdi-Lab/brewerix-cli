from itertools import groupby
from operator import itemgetter
from typing import List, Iterator, Tuple

VERSION = 1


# logfile = open("memory-profiler.log", "w+")


# @profile(stream=logfile)
def create_guess_loi_table(lines: Iterator[list], head: List, output: str = "final-output-table.txt"):
    # reduced_snps = sort_by_columns(reduce_snp_redundancies(lines, gene_col), [0, 5, 1])  # sort by chr, gene, position
    reduced_snps = sort_by_columns(lines, [0, 5, 1])  # sort by chr, gene, position
    write_guess_loi_table(reduced_snps, head, output)


def sort_by_columns(lines: Iterator[List], columns: List) -> List:
    yield from sorted(lines, key=itemgetter(*columns))


def reduce_snp_redundancies(snp_lines: Iterator, gene_col: int) -> List:
    # DEPRECATED
    for key, values in groupby(snp_lines, itemgetter(gene_col)):
        # TODO: use the iterator rather than the list?
        values = list(values)
        interesting_snps, overall_gene = extract_informative_snps(values, gene_col, 0.1)

        if len(interesting_snps) == 0:
            first_line_gene_annotation = values[0][:gene_col + 1]
            yield collapse_to_gene_info(first_line_gene_annotation, overall_gene)
        else:
            yield from interesting_snps


def extract_informative_snps(values: List, gene_col: int, ratio_min: float = 0.1) -> Tuple:
    # DEPRECATED

    interesting_snps = []
    overall_gene_expression = []

    for line in values:
        allelic_expressions = line[(gene_col + 1):]
        a_ratio = [ratio_allele_expression(x) for x in allelic_expressions]
        a_sum = [sum_allele_expression(x) for x in allelic_expressions]
        overall_gene_expression.append(a_sum)

        for r in a_ratio:
            if r >= ratio_min:
                interesting_snps.append(line)
                break

    overall_gene = [sum(x) for x in zip(*overall_gene_expression)]
    return interesting_snps, overall_gene


def collapse_to_gene_info(gene_annotation: List, overall_gene_expression: List, snp_id_text: str = "rs_multi",
                          fake_pvalue: float = 1.0) -> List:
    gene_annotation[2] = snp_id_text
    # gene_annotation = [str(x) for x in gene_annotation]

    ref_alt_fake = zip(overall_gene_expression, [0] * len(overall_gene_expression))
    ref_alt_fake = add_fake_pvalue(ref_alt_fake, value=fake_pvalue)

    overall_gene_expression = [','.join([str(v) for v in x]) for x in ref_alt_fake]
    return gene_annotation + overall_gene_expression


def add_fake_pvalue(l: Iterator, value=1.0):
    for t in l:
        yield [v for v in t] + [value]


def write_guess_loi_table(lines: Iterator[List], header: List, filename: str) -> None:
    with open(filename, 'wt') as out:
        out.write('# version=%d\n' % VERSION)
        out.write('\t'.join(header) + '\n')
        for line in lines:
            line_str = [str(x) for x in line]
            out.write('\t'.join(line_str) + '\n')


def sum_allele_expression(ae):
    if ae == "NA":
        return 0

    aes = ae.split(',')[:2]  # exclude the pvalues
    return sum([int(x) for x in aes])


def ratio_allele_expression(ae):
    if ae == 'NA':
        return 0

    aes = ae.split(',')[:2]  # exclude the pvalues
    aes = [int(x) for x in aes]
    return min(aes) / float(sum(aes))
