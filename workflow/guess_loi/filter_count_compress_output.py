from collections import namedtuple
from itertools import islice
from typing import Tuple, List

from workflow.guess_loi.table import sort_by_columns, sum_allele_expression, ratio_allele_expression

GeneSnpSummary = namedtuple('GeneSnpSummary', 'snps, gene')


def sort_file_by_gene_name_and_position(lines) -> list:
    # header, lines = read_annoted_ase(in_file)
    # header = next(lines)
    # return header, sort_by_columns(lines, [5, 1])
    return sort_by_columns(lines, [5, 1])


def read_annoted_ase(in_file: str) -> Tuple[str, List]:
    # DEPRECATED
    lines = []
    with open(in_file, "rt") as fd:
        header = next(fd)

        for idx, line in enumerate(fd, 2):
            tks = line.rstrip('\n').split('\t')
            lines.append(tks)
        return header, lines


def increment_values(gene, annotation, values, goe):
    if gene in goe:
        goe[gene][2] = [sum(x) for x in zip(goe[gene][2], values)]
    else:
        goe[gene] = [annotation, [gene], values]


def compute_overall_expression(snp_lines, gene_col):
    gene_overall_expression = {}
    for line in list(islice(snp_lines, 1000)):
        allelic_expressions = line[(gene_col + 1):]
        gene_name = line[gene_col]
        annotation = line[:gene_col]

        a_sum = [sum_allele_expression(x) for x in allelic_expressions]
        increment_values(gene_name, annotation, a_sum, gene_overall_expression)

    return gene_overall_expression


def filter_useful_snps(snp_lines, gene_col, ratio_min):
    for line in snp_lines:
        allelic_expressions = line[gene_col+1:]
        a_ratio = [ratio_allele_expression(x) for x in allelic_expressions]
        for r in a_ratio:
            if r >= ratio_min:
                yield(line)
                break
