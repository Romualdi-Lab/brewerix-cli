from collections import namedtuple
from itertools import groupby
from operator import itemgetter
from typing import Tuple, List, Iterator

GeneSnpSummary = namedtuple('GeneSnpSummary', 'snps, gene')


def reduce_snp_redundancies(snp_lines: Iterator, gene_col: int) -> List:
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


def collapse_to_gene_info(gene_annotation: List, overall_gene_expression: List, snp_id_text: str="rs_multi") -> List:
    gene_annotation[2] = snp_id_text
    ref_alt_fake = zip(overall_gene_expression, [0 for _ in range(len(overall_gene_expression))])
    overall_gene_expression = [','.join([str(v) for v in x]) for x in ref_alt_fake]
    return gene_annotation + overall_gene_expression


def sort_file_by_gene_name_and_position(lines) -> Tuple[str, list]:
    # header, lines = read_annoted_ase(in_file)
    header = next(lines)
    return header, sort_by_columns(lines, [5, 1])


def read_annoted_ase(in_file: str) -> Tuple[str, List]:
    lines = []
    with open(in_file, "rt") as fd:
        header = next(fd)

        for idx, line in enumerate(fd, 2):
            tks = line.rstrip('\n').split('\t')
            lines.append(tks)
        return header, lines


def sort_by_columns(lines: List, columns: List) -> List:
    yield from sorted(lines, key=itemgetter(*columns))


def write_guess_loi_table(lines: Iterator[List], header: List, filename: str) -> None:
    with open(filename, 'wt') as out:
        out.write('\t'.join(header) + '\n')
        for line in lines:
            out.write('\t'.join(line) + '\n')


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
