from tempfile import NamedTemporaryFile

from intervaltree import IntervalTree
from pytest import raises

from workflow.guess_loi.snp_gene_association import AnnotationError, read_bed_index, annotate_aser, read_ase_table


def test_fail_read_bed_unsorted():
    with NamedTemporaryFile("w+t") as f:
        content = [
            "1\t2\t10\tgene_a\t+",
            "2\t6\t50\tgene_b\t+",
            "1\t30\t45\tgene_c\t-",
        ]
        f.write('\n'.join(content))
        f.flush()

        with raises(AnnotationError):
            list(read_bed_index(f.name))


def test_fail_read_bed_index_with_less_than_5_fileds():
    with NamedTemporaryFile("w+t") as f:
        content = [
            "1\t2\t10\tgene_a",
            "2\t6\t50\tgene_b",
            "1\t30\t45\tgene_c",
        ]
        f.write('\n'.join(content))
        f.flush()

        with raises(AnnotationError):
            list(read_bed_index(f.name))


def test_fail_read_ase_table_malformed_id():
    with NamedTemporaryFile("w+t") as f:
        content = [
            "chr_pos_rs_ref_alt\ts1\ts2",
            "1_20_rs123\t1/2\t2/2",
        ]
        f.write('\n'.join(content))
        f.flush()

        with raises(AnnotationError):
            list(read_ase_table(f.name))


def test_annotate_aser_simple():
    infos_head = ['chr', 'pos', 'rs', 'ref', 'alt']
    values_head = ['sample1', 'sample2']
    
    infos = ['1', '20', 'rs123', 'A', 'T']
    values = ['r1/a1', 'r2/a2']
    
    lines = [
        ['chr', 'pos', infos_head, values_head],
        ['1', '20', infos, values],
    ]
    tree = IntervalTree()
    tree[5:300] = "Gene_a"
    bed_idx = {'1': tree}

    annotated = [
        infos_head + ['symbol'] + values_head,
        infos + ['Gene_a'] + values,
    ]

    assert list(annotate_aser(lines, bed_idx)) == annotated
