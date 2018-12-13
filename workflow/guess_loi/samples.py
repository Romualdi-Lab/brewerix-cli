import re
from collections import namedtuple, defaultdict
from os.path import basename


Sample = namedtuple('Sample', 'name, fastqs, bam')


def single_samples(fastqs):
    rx = re.compile(r'(.*)\.f(ast)?q(\.gz)?$')

    samples = []
    for fastq in fastqs:
        m = rx.match(fastq)
        if m is None:
            exit('invalid filename for single-end reads: ' + fastq)
        else:
            samples.append(Sample(basename(m.group(1)), [fastq], None))

    return samples


def paired_samples(fastqs):
    rx = re.compile(r'(.*)_[12]\.f(ast)?q(\.gz)?$')

    groups = defaultdict(list)
    for fastq in fastqs:
        m = rx.match(fastq)
        if m is None:
            exit('invalid filename for paired-end reads: ' + fastq)
        else:
            groups[basename(m.group(1))].append(fastq)

    samples = []
    for name, files in groups.items():
        if len(files) != 2:
            exit('wrong number of files for group %r: %d' % (name, len(files)))
        else:
            samples.append(Sample(name, sorted(files), None))

    return samples
