import re
from itertools import groupby
from operator import itemgetter
from sys import argv
from typing import Dict, List, TextIO, Iterable, Iterator


class AnnotationError(Exception):
    pass


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
    bedidx = read_bed_index(bed_file)

    with open(vcf_file, 'r') as fd:
        for idx, line in enumerate(fd):
            if line[0] == '#':
                continue
            chromosome, pos, rs_id, ref, alt, *infos = line.rstrip('\n').split('\t')
            snp_internal_id = '_'.join([chromosome, pos, rs_id, ref, alt])

            if chromosome in bedidx:
                for start, stop, gene, _strand in bedidx[chromosome]:
                    if start <= pos <= stop:

                        print(snp_internal_id + '\t' + gene)


def annotate():
    aser_table = argv[1]
    bed_file = argv[2]
    annotate_aser_table_from_bed(aser_table, bed_file)


def annotate_aser_table_from_bed(aser_table_file: str, bed_file: str) -> Iterator:
    # the first column is the id. The rest are the values.
    bed_idx = read_bed_index(bed_file)
    # TODO: remove pseudoautosomal region
    read_lines = read_ase_table(aser_table_file)
    annotated_lines = annotate_aser(read_lines, bed_idx)
    return annotated_lines
    # write_annotated_aser_table(annotated_lines, output)


def read_bed_index(bed_file: str) -> Dict:
    with open(bed_file, 'rt') as fd:
        return read_bed(fd)


def read_bed(fd: TextIO) -> Dict:
    return {chromosome: [line[1:] for line in grp] for chromosome, grp in groupby(read_line_bed(fd), itemgetter(0))}


def read_line_bed(fd: TextIO) -> Iterable[List]:
    # TODO: check if bed has 5 tokens
    pre_id = None

    for lineno, line in enumerate(fd, 1):
        tokens = line.rstrip().split('\t')

        if len(tokens) != 5:
            raise AnnotationError("malformed input: incorrect number of columns at line %s" % lineno)

        if pre_id is not None and tokens[0] < pre_id:
            raise AnnotationError("malformed input: lexicographically sorted on col 1 at line %s" % lineno)

        pre_id = tokens[0]

        yield tokens


def read_ase_table(filename: str) -> Iterable[List]:
    with open(filename, 'r') as fd:
        length_t = None
        for lineno, line in enumerate(fd):
            tks = line.rstrip('\n').split('\t')
            l_err = lineno + 1
            if length_t is not None and length_t != len(tks):
                raise AnnotationError("malformed input: incorrect number of columns at line %d" % l_err)
            length_t = len(tks)
            infos = tks[0].split("_")

            if len(infos) != 5:
                raise AnnotationError("malformed input: id must be format as %r at line %d" %
                                      ('chr_pos_rs_ref_alt', l_err))

            yield lineno, infos[0], infos[1], infos, tks[1:]


def annotate_aser(lines: Iterable[List], bed_idx: Dict) -> List:
    for lineno, chrom, pos, infos, values in lines:
        if lineno == 0:
            yield infos + ["symbol"] + values

        if chrom in bed_idx:
            for start, stop, gene, _strand in bed_idx[chrom]:
                if int(start) <= int(pos) <= int(stop):
                    yield infos + [gene] + values


def write_annotated_aser_table(lines: Iterable[str], filename):
    with open(filename, 'w+t') as f:
        for line in lines:
            f.write(line + '\n')


if __name__ == '__main__':
    annotate()
