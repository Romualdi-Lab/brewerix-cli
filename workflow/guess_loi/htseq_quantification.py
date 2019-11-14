from argparse import ArgumentParser
from subprocess import check_call

# htseq-count -s no -f bam HPD01_naive.bam human-reduced.gtf >HPD01_naive-no-str.quant.txt


def quantify():
    parser = ArgumentParser(description="""
                Quantify bam files using htseq-count
                """)
    parser.add_argument('gtf', help="gtf file name")
    parser.add_argument('bams', nargs='+', help="bam files")
    parser.add_argument('--no-head', dest='nohead', help="increase output verbosity",
                        action="store_false", default=True)
    parser.add_argument('-s', '--strandness', dest='strand', help="either no, yes or reverse (TrueSeq)", default="no")
    parser.add_argument('-r', '--order', dest='order', help="either pos or name; default pos", default="pos")

    args = parser.parse_args()

    run_quantification(args.bams, args.gtf, args.nohead, args.strand, args.order)


def run_quantification(files, gtf, head, strand, order):
    header = ["gene"] + files
    if head:
        print('\t'.join(header))

    cmd = ['htseq-count',
           '-f', 'bam',
           "-r", order,
           '-s', strand] + files + [gtf]

    check_call(cmd)
