from setuptools import find_packages, setup

setup(
    name="workflow-guess-loi",
    version="0.1.0",
    author="Romualdi's Lab",
    author_email=[
        "paolo.cavei@gmail.com",
        "enrica.calura@gmail.com",
        "gbrsales@gmail.com",
    ],
    description="Pipeline for Loss Of Imprinting Identification",
    long_description="""
            Collection of tools to perform loss of imprinting analysis.
        """,
    url="https://gitlab.romualdi.bio.unipd.it:workflow/guess_loi.git",
    packages=find_packages(),
    entry_points={
        'console_scripts': [
            'guess_loi=workflow.guess_loi.guess_loi:guess_loi_from_fqs',
            'format_ase=workflow.guess_loi.format_ase:format_ase',
            'annotate=workflow.guess_loi.snp_gene_association:annotate',
            'compact-snp=workflow.guess_loi.filter_count_compress_output:compact_snps',
            'hisat2-buildIndex=workflow.guess_loi.hisat_build_index:hisat_build_index',
            'run-hisat2=workflow.guess_loi.hisat2_align:run_hisat2',
        ],
    }
)


