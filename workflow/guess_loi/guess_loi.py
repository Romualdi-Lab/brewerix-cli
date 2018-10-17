from os.path import expanduser
from subprocess import call, DEVNULL
from sys import argv

from workflow.guess_loi.format_ase import format_ase_internal


def guess_loi():
    print(argv)
    bams = argv[3:]
    SNPs = argv[1]
    bed = argv[2]
    samples = [ bam.rstrip(".bam") for bam in bams ]

    print(bams)
    print(SNPs)
    print(bed)
    print(samples)

    gatk = check_gatk(gatk="~/local/stow/gatk-4.0.4.0/gatk")

    # check Genome Index
    # sort bam if need be
    # create RG tags, if need be

    aser_count(gatk, bams, SNPs, bed, samples)

    for sample in samples:
        intermediate_file = ''.join([sample, ".ASER.txt"])
        with open(intermediate_file, 'r') as fd:
            intermediate_out =''.join([sample, ".aser"])
            with open(intermediate_out, "w") as intermediate:
                format_ase_internal(fd, intermediate)


def check_gatk(gatk='~/local/bin/gatk'):
    try:
        gatk = gatk.replace("~", expanduser("~"))
        call([gatk, "--list"],
             stdin=DEVNULL, stdout=DEVNULL, stderr=DEVNULL)
    except FileNotFoundError as e:
        print("check if gatk is installed.")
        exit(134)

    return(gatk)

def aser_count(gatk, bams, SNPs, bed, samples):
    # run aser count
    # format the results
    # annotate with bed file
    # create allele count info (compress the output)
    for bam, sample in zip(bams, samples):
        run_aser_on_bam(gatk, bam, SNPs, "genome-sort.fa", sample)



def run_aser_on_bam(gatk, bam, vcf, genome, sample):
    call([gatk, "ASEReadCounter",
          "-I", bam,
          "-V", vcf,
          "-R", genome,
          "-O", ''.join([ sample,".ASER.txt"])])


if __name__ == '__main__':
    guess_loi()

