from tempfile import NamedTemporaryFile

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
    lines = [
        [0, 'chr', 'pos', 'chr\tpos\trs\tref\talt', 'sample1\tsample2'],
        [1, '1', '20', '1\t20\trs123\tA\tT', 'r1/a1\tr2/a2'],
    ]
    bed_idx = {'1': [[5, 300, "Gene_a", '+']]}

    annotated = [
        'chr\tpos\trs\tref\talt\tsymbol\tsample1\tsample2',
        '1\t20\trs123\tA\tT\tGene_a\tr1/a1\tr2/a2'
    ]

    assert list(annotate_aser(lines, bed_idx)) == annotated
