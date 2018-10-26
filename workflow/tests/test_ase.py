from tempfile import NamedTemporaryFile
from typing import List

from pytest import raises

from workflow.guess_loi.ase import AseError, Entry, format_ase, read_gatk_ase, merge_ase, Ase


def test_failure_with_non_int_values():
    with raises(AseError):
        entry = [Entry("id", "a", 0, 1)]
        list(format_ase(entry, 5, "filein"))


def test_failure_file_with_no_header():
    with NamedTemporaryFile('w+t') as f:
        f.write("There is no header\n")
        f.flush()

        with raises(AseError):
            list(read_gatk_ase(f.name))

def test_merge():
    ases = [
        {
            'k1': 'r1/a1',
            'k2': 'r2/a2',
        },
        {
            'k1': 'r3/a3',
            'k3': 'r4/a4',
        },
    ]  # type: List[Ase]

    merged = {
        'k1': 'r1/a1\tr3/a3',
        'k2': 'r2/a2\tNA',
        'k3': 'NA\tr4/a4',
    }

    assert merge_ase(ases) == merged
