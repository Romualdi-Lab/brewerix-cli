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
            'guess_loi=workflow.guess_loi.guess_loi:guess_loi',
            'format_ase=workflow.guess_loi.format_ase:format_ase',
            'snp2gene=workflow.guess_loi.snp_gene_association:snp2gene',
        ],
    }
)

