from os.path import abspath, join
from subprocess import check_call
from tempfile import TemporaryDirectory

import pytest

from workflow.guess_loi.alignments import align
from workflow.guess_loi.samples import Sample
from workflow.guess_loi.samtools import call_samtools_index

genome_fa = abspath("./workflow/testfiles/genome.fa")
fastq_file = abspath("./workflow/testfiles/single_trueSeq.fq.gz")
index = abspath("./workflow/testfiles/index/index")
bed = abspath("./workflow/testfiles/genes_regions.bed")

@pytest.fixture(scope="session")
def hisat2_genome_index():
    with TemporaryDirectory() as wdir:
        check_call(['hisat2-build', '--quiet', genome_fa, 'index'], cwd=wdir)
        yield join(wdir, 'index')


@pytest.fixture(scope="session")
def hisat2_bam():
    with TemporaryDirectory() as wdir:
        s = Sample('single_trueSeq', [fastq_file], "single_trueSeq.bam")
        align(s, index, bed, s.bam, 1, 1)
        call_samtools_index(s.bam)
        yield join(wdir, s.bam)
