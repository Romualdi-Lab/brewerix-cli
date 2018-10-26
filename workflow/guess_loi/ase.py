from collections import namedtuple
from subprocess import check_call
from typing import Tuple, Iterable, Dict, List, NewType

from workflow.guess_loi.filter_count_compress_output import reduce_snp_redundancies, \
    write_guess_loi_table, sort_by_columns

Entry = namedtuple('Entry', 'gene_id, ref, alt, lineno')

AseKey = NewType('AseKey', str)

Ase = Dict[AseKey, str]

Sample = NewType('Sample', str)


class AseError(Exception):
    pass


def ase_table(gatk, bams, snps, genome, samples: List[Sample]) -> str:
    output = "ASER_table.txt"
    ases = [aser_count(gatk, bam, snps, genome, sample) for bam, sample in zip(bams, samples)]
    write_ase(merge_ase(ases), samples, output)
    return output


def aser_count(gatk, bam, vcf, genome, sample):
    step1 = sample + ".ASER.txt"

    check_call([
        gatk, "ASEReadCounter",
        "-I", bam,
        "-V", vcf,
        "-R", genome,
        "-O", step1,
    ])

    return read_ase(step1)


def read_ase(input_: str, min_coverage: int = 5) -> Ase:
    entries = read_gatk_ase(input_)
    return dict(format_ase(entries, min_coverage, input_))


def read_gatk_ase(input_: str) -> Iterable[Entry]:
    with open(input_, 'rt') as fd:
        line = next(fd)
        if not line.startswith('contig'):
            raise AseError('missing header in file %r' % input_)

        for lineno, line in enumerate(fd, 2):
            tokens = line.rstrip('\n').split('\t')
            if len(tokens) < 7:
                raise AseError("less than 7 tokens at line %d in file %r" % (lineno, input_))

            yield Entry('_'.join(tokens[:5]), tokens[5], tokens[6], lineno)


def format_ase(entries: Iterable[Entry], min_coverage: int, filename: str) -> Iterable[Tuple[AseKey, str]]:
    for entry in entries:
        try:
            ref_i = int(entry.ref)
            alt_i = int(entry.alt)
        except ValueError:
            raise AseError('invalid ref/alt at line %d in file %r' % (-1, filename))

        if ref_i + alt_i >= min_coverage:
            ref_alt = entry.ref + ',' + entry.alt
            yield entry.gene_id, ref_alt


def merge_ase(ases: List[Ase]) -> Ase:
    all_keys = merge_keys(ases)
    return gather_values(all_keys, ases)


def merge_keys(ases: List[Ase]) -> List[AseKey]:
    keys = set()
    for ase in ases:
        keys |= set(ase.keys())
    return sorted(keys)


def gather_values(all_keys: List[AseKey], ases: List[Ase]) -> Ase:
    return {key: '\t'.join([ase.get(key, "NA") for ase in ases]) for key in all_keys}


def write_ase(ase: Ase, samples: List[Sample], output: str) -> None:
    with open(output, 'wt') as fd:
        fd.write("chr_pos_rs_ref_alt\t" + '\t'.join(samples) + '\n')

        for key, value in ase.items():
            fd.write('%s\t%s\n' % (key, value))


def create_guess_loi_table(lines: List, head: List, gene_col: int=5, output: str= "final-output-table.txt"):
    reduced_snps = sort_by_columns(reduce_snp_redundancies(lines, gene_col), [0, 1])
    write_guess_loi_table(reduced_snps, head, output)


def create_keys_dictionaries(files):
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
