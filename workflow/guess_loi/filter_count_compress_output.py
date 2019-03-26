from collections import namedtuple
from typing import Tuple, List

from workflow.guess_loi.table import sort_by_columns

GeneSnpSummary = namedtuple('GeneSnpSummary', 'snps, gene')


def sort_file_by_gene_name_and_position(lines) -> Tuple[List, list]:
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


