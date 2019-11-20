from os.path import join
from re import search
from subprocess import PIPE, check_output, check_call, Popen
from tempfile import TemporaryDirectory


def call_samtools_index(bam):
    check_call(["samtools", "index", bam])


def check_rg_tag(bam_input):
    out = check_output(["samtools", "view", "-H", bam_input], universal_newlines=True)
    for line in out.split("\n"):
        if not search(r'^@RG', line):
            break
        return True

    return False


def samtools_view(stdin, bed):
    args = [
        'samtools',
        "view",
        "-u",
        "-L",
        bed
    ]

    return Popen(args, stdin=stdin, stdout=PIPE)


def split_bam_by_chromosomes(bam, wdir, chromosomes):
    filenames = []
    for chrom in chromosomes:
        filename = join(wdir, chrom + '_' + bam)
        filenames.append(filename)
        check_call([
            'samtools',
            "view",
            '-u',
            bam,
            chrom,
            '-o', filename
        ])
        call_samtools_index(filename)
    return filenames


def split_vcf_by_chromosomes(vcf, wdir, chromosomes):
    filenames = []
    with TemporaryDirectory() as twdir:
        vcf_gz = join(twdir, vcf + '.gz')
        with open(vcf_gz, "wb") as fd:
            check_call(["bgzip", "-c", vcf], stdout=fd)

        check_call(["tabix", "-p", "vcf", vcf_gz])

        for chrom in chromosomes:
            filename = join(wdir, chrom + '_' + vcf)
            filenames.append(filename)
            check_call([
                'bcftools',
                "view",
                vcf_gz,
                chrom,
                '-o', filename
            ])
            check_call('gatk IndexFeatureFile -F ' + filename, shell=True)

        return filenames


def samtools_sort(stdin, output, samtools_threads):
    args = [
        'samtools',
        'sort',
        '-l', '9',
        '-m', '1G',
        '-@', str(samtools_threads),
        '-o', output
    ]

    return Popen(args, stdin=stdin)
