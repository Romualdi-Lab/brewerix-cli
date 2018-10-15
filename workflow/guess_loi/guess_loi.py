from subprocess import call


def guess_loi():
    bams = ["1.bam", "2.bam"]
    # check_gatk()
    # check Genome Index
    # sort bam if need be
    # create RG tags, if need be

    aser_count(bams, SNPs, bed)

    print(bams)


def check_gatk(gatk='~/local/bin/gatk'):
    try:
        call([gatk, "--list"])
    except FileNotFoundError as e:
        print("check if gatk is installed.")
        # exit(134)

def aser_count(bams, SNPs, bed, sample):
    # run aser count
    # format the results
    # annotate with bed file
    # create allele count info (compress the output)
    for bam in bams:
        run_aser_on_bam(bam, SNPs, bed, sample)



def run_aser_on_bam(bam, vcf, genome, sample):
    call(["gatk", "ASEReadCounter",
          "-I", bam,
          "-V", vcf,
          "-R", genome,
          "-O", ''.join(sample,".ASER.txt")])


if __name__ == '__main__':
    guess_loi()

