import argparse


def parse_args():
    parser = argparse.ArgumentParser(description="""
    Create Tables for BrewerIX
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
    parser_bams.add_argument('--gatk-memory', dest="gatkmem", type=int, default=None,
                             help="Integer, the number of Gb to give to GATK")
    parser_bams.add_argument('--remove-duplicates', dest="rm_dup", default=False, action='store_true',
                            help="Use Picard Tools to mark duplicates. Mark in the very same bam file.")
    parser_bams.add_argument('--clean', dest="clean", default=False, action='store_true',
                            help="clean intermediate files")

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
    parser_fqs.add_argument('--gatk-memory', dest="gatkmem", type=int, default=None,
                             help="Integer, the number of Gb to give to GATK")
    parser_fqs.add_argument('--remove-duplicates', dest="rm_dup", default=False, action='store_true',
                            help="Use Picard Tools to mark duplicates")
    parser_fqs.add_argument('--split-N-sigars', dest="split_N", default=False, action='store_true',
                            help="Use Gatk Tools to split N sigars.")
    parser_fqs.add_argument('--do-recalibration', dest="do_recal", default=False, action='store_true',
                            help="Use Gatk Tools to perform recalibration. This activate both split-N-sigar "
                                 "and remove duplicates options")
    parser_fqs.add_argument('--clean', dest="clean", default=False, action='store_true',
                            help="clean intermediate files")

    return parser.parse_args()
