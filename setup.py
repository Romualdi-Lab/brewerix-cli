from setuptools import find_packages, setup

setup(
    name="brewerix-cli",
    version="0.9.0",
    author="Romualdi's Lab",
    author_email=[
        "paolo.cavei@gmail.com",
        "gbrsales@gmail.com",
    ],
    description="Pipeline for Loss Of Imprinting and X Chromosome Inactivation Analysis",
    long_description="""
            Collection of tools to perform loss of imprinting and X chromosome inactivation analysis.
        """,
    url="https://gitlab.romualdi.bio.unipd.it:workflow/guess_loi.git",
    packages=find_packages(),
    install_requires=[
        'intervaltree',
        'scipy',
        'pysam'
    ],
    extras_require={
        'test': [
            'pytest',
        ],
    },
    entry_points={
        'console_scripts': [
            'guess-LOI=workflow.guess_loi.guess_loi:guess_loi',
            # 'haplotype_caller=workflow.guess_loi.haplotype_caller_rna:haplotype_caller',
            # 'haplotype_caller_dna=workflow.guess_loi.haplotype_caller_dna:haplotype_caller_dna',
            # 'select_variants=workflow.guess_loi.select_variants:select_variants',
            # 'merge_vcfs=workflow.guess_loi.merge_vcfs:merge_vcfs',
            # 'concat_vcfs=workflow.guess_loi.concat_vcfs:concat_vcfs',
            # 'call_and_merge=workflow.guess_loi.add_haplo_callad_from_multiallelic:call_and_merge',
            # 'het_genotype=workflow.guess_loi.vcf_related_functions:reformat_vcf_with_het_genotype',
            # 'annotate=workflow.guess_loi.snp_gene_association:annotate',
            # 'hisat2-buildIndex=workflow.guess_loi.hisat_build_index:hisat_build_index',
            # 'run-hisat2=workflow.guess_loi.hisat2_align:run_hisat2',
            # 'fake_genotype=workflow.guess_loi.vcf_related_functions:reformat_vcf_with_het_genotype',
            # 'exclude_invalid_snv_ids=workflow.guess_loi.vcf_related_functions:exclude_invalid_snv_ids',
        ],
    }
)
