from os.path import exists, basename, join
from sys import argv
from tempfile import TemporaryDirectory

from workflow.guess_loi.checks import check_command_availability
from workflow.guess_loi.haplotype_caller_rna import run_haplotype_caller
from workflow.guess_loi.merge_vcfs import run_merge_vcfs
from workflow.guess_loi.select_variants import run_select_variants


def call_and_merge():
    check_command_availability(['gatk'])
    if len(argv) == 1 or len(argv) < 8:
        print("Run like this:")
        print('%s OUTPUT out_haplo_vcf out_selector_vcf biallelic_vcf multiallelic_vcf genome bams...\n' % basename(argv[0]))
        exit(0)

    output_vcf = argv[1]
    output_haplo_caller = argv[2]
    output_selector = argv[3]
    biallelic_snp = argv[4]
    multiallelic_snp = argv[5]
    genome = argv[6]
    bams = argv[7:]

    call_and_merge_multiallelic_to_biallelic(output_vcf, output_haplo_caller, output_selector, biallelic_snp,
                                             multiallelic_snp, bams, genome)


def call_and_merge_multiallelic_to_biallelic(output_vcf, output_haplo_caller, output_selector, biallelic_snp,
                                             multiallelic_snp, bams, genome):

    if not exists(output_haplo_caller):
        print("Running Haplotype Caller...")
        run_haplotype_caller(multiallelic_snp, bams, genome, output_haplo_caller)

    if not exists(output_selector):
        print("Running Selector...")
        run_select_variants(genome, output_haplo_caller,
                            ["--select-type-to-exclude", "INDEL", "--restrict-alleles-to", "BIALLELIC"],
                            output_selector)

    print("Merging files...")
    with TemporaryDirectory() as wdir:
        intermediate_file = join(wdir, output_vcf)
        run_merge_vcfs([biallelic_snp, output_selector], intermediate_file)
        run_select_variants(genome, intermediate_file,
                            ["--select-type-to-exclude", "INDEL", "--restrict-alleles-to", "BIALLELIC"], output_vcf)
