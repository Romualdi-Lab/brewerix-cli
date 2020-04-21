from argparse import ArgumentParser
from subprocess import check_call

from workflow.guess_loi.checks import check_command_availability


def feature_counts():
    parser = ArgumentParser(description="""
                Quantify bam files using featureCount from 'http://subread.sourceforge.net/'
                """)
    parser.add_argument('outfile', help="output file name")
    parser.add_argument('gtf', help="gtf file name")
    parser.add_argument('bams', nargs='+', help="bam files")
    parser.add_argument('-s', '--strandness', dest='strand',
                        help="either 0, 1, 2 for unstranded, stranded, reverse (TrueSeq) "
                             "respectively. A line of comma separated for each sample",
                        default="0")
    parser.add_argument('-p', '--paired', action='store_true', help="for paired end", default=False)
    parser.add_argument('-t', '--threads', type=int, help="number of cpus", default=1)
    parser.add_argument('-d', '--dry-run', action="store_true", help="do not run the command", default=False)

    args = parser.parse_args()
    # featureCounts -T 32 -p -s 2 -a human-ens94.gtf -o raw-counts-fc-reverse.txt *.bam
    run_quantification(args.bams, args.gtf, args.output, args.strand, args.paired, args.threads, args.dry_run)


def run_quantification(files, gtf, output, strand, paired, threads, dry_run):
    check_command_availability(['featureCounts'])
    strandness = {'0': "unstranded",
                  '1': "stranded",
                  '2': "reverse"}

    paired = '-p' if paired else ''
    cmd = ['featureCounts',
           '-T', threads,
           paired,
           '-s', strand,
           '-F', "gtf",
           '-a', gtf,
           '-o', output,
           ] + files
    print("Run feautureCounts in %s" % strandness[strand])
    print(' '.join(cmd))
    if not dry_run:
        check_call(cmd)
