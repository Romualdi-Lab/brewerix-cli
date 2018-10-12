from subprocess import call


def guess_loi():
    bams = ["1.bam", "2.bam"]
    check_gatk()

    print(bams)

def check_gatk(gatk='~/local/bin/gatk'):
    try:
        call([gatk, "--list"])
    except FileNotFoundError as e:
        print("check if gatk is installed.")
        # exit(134)



if __name__ == '__main__':
    guess_loi()

