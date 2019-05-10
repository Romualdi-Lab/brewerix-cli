from workflow.guess_loi.guess_loi import create_annotated_lines


def test_filter_useful_snps():
    goe = {'A': [["chr1", 1, "rs_A"], ["A"], [3, 6, 180]],
           'B': [["chr1", 4, "rs_B"], ["B"], [1, 1, 10]],
           'C': [["chr1", 4, "rs_C"], ["C"], [2, 2, 10]],
           'D': [["chr2", 4, "rs_D"], ["D"], [5, 6, 10]]}

    good_snps = [["chr1", 1, "rs_A", "A", '4,5,0.05', '2,0,1', '60,0,1.0'],
                 ["chr1", 4, "rs_B", "B", '1,0,1.0', '21,20,0.001', '10,0,1'],
                 ["chr2", 9, "rs_C", "C", '115,100,0.001', '2,0,1', '60,0,1']]

    genes2tss = {"A": 1, "B": 3, "C": 4, "D": 2}

    good_out = []
    for line in good_snps:
        out = line[:]
        out.insert(4, int(genes2tss[line[3]]))
        good_out.append(out)

    good_out.append(['chr2', 4, 'rs_multi', 'D', 2, '5,0,1.0', '6,0,1.0', '10,0,1.0'])

    assert list(create_annotated_lines(good_snps, goe, gene_col=3, genes2tss=genes2tss)) == good_out
