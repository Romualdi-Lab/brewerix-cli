from workflow.guess_loi.filter_count_compress_output import increment_values, compute_overall_expression, \
    filter_useful_snps


def test_increment_values():
    goe = {}
    gene = "A"
    annotation = ["ch1", "1456"]
    values = [4, 5]

    increment_values(gene, annotation, values, goe)
    increment_values(gene, annotation, values, goe)

    final_out = {'A': [annotation, [gene], [8, 10]]}

    assert goe == final_out


def test_compute_overall_expression():
    snp_lines = [["chr1", 1, "A", '1,0,1', '2,0,1', '60,0,1'],
                 ["chr1", 2, "A", '1,0,1', '2,0,1', '60,0,1'],
                 ["chr1", 4, "B", '1,0,1', '1,0,1', '10,0,1'],
                 ["chr1", 3, "A", '1,0,1', '2,0,1', '60,0,1']]

    gene_col = 2

    goe = {'A': [["chr1", 1], ["A"], [3, 6, 180]],
           'B': [["chr1", 4], ["B"], [1, 1, 10]]
           }

    assert compute_overall_expression(snp_lines, gene_col) == goe


def test_filter_useful_snps():
    snp_lines = [["chr1", 1, "A", '4,5,0.05', '2,0,1', '60,0,1'],
                 ["chr1", 1, "A", '4,0,1', '2,0,1', '1,11,0.9'],
                 ["chr1", 4, "B", '1,0,1', '21,20,0.001', '10,0,1'],
                 ["chr2", 8, "C", '1,0,1', '2,0,1', '60,0,1'],
                 ["chr2", 9, "C", '115,100,0.001', '2,0,1', '60,0,1']]

    good_snps = [["chr1", 1, "A", '4,5,0.05', '2,0,1', '60,0,1'],
                 ["chr1", 4, "B", '1,0,1', '21,20,0.001', '10,0,1'],
                 ["chr2", 9, "C", '115,100,0.001', '2,0,1', '60,0,1']]

    assert list(filter_useful_snps(snp_lines, gene_col=2, ratio_min=0.1)) == good_snps
