import argparse


def parse_args():
    parser = argparse.ArgumentParser(description="""
    Create Tables for GUESS LOY Tool
    """)

    subparsers = parser.add_subparsers(dest='mode')

    parser_bams = subparsers.add_parser('bams', help='stats from bams')

    parser_bams.add_argument('snps', help="vcfs of the SNPs. Index must be there as well")
    parser_bams.add_argument('bed', help="the bed with the boundaries of the imprinted genes")
    parser_bams.add_argument('genome_dict', help="the genome dict. The reference fa")
    parser_bams.add_argument('bams', nargs='+', help="the bams")
    parser_bams.add_argument('--threads', type=int, default=1,
                             help="Number of CPUs to use; please note that any CPU use approximately 2G RAM")
    parser_bams.add_argument('--write-progress-to', dest="progress",
                             help="write analysis progress to file")
    parser_bams.add_argument('--multi', help="the vcf files with multi alleles")

    parser_fqs = subparsers.add_parser('fqs', help='stats from fastqs')

    parser_fqs.add_argument('snps', help="vcfs of the SNPs. Index must be there as well")
    parser_fqs.add_argument('bed', help="the bed with the boundaries of the imprinted genes")
    parser_fqs.add_argument('genome_dict', help="the genome dict. The reference fa and dict")
    parser_fqs.add_argument('genome_idx', help="the genome index for hisat2")
    parser_fqs.add_argument('fqs', nargs='+', help="the fqs")
    parser_fqs.add_argument('-p', '--paired', dest='is_paired', default=False,
                            action='store_true', help="paired end alignment")
    parser_fqs.add_argument('--threads', type=int, default=1,
                            help="Number of CPUs to use; please note that any CPU use approximately 2G RAM")
    parser_fqs.add_argument('--write-progress-to', dest="progress",
                            help="write analysis progress to file")
    parser_fqs.add_argument('--multi', help="the vcf files with multi alleles")

    return parser.parse_args()
