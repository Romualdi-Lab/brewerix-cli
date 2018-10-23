from os.path import exists, expanduser
from subprocess import call, DEVNULL
from tempfile import TemporaryFile


def check_file_exists(file):
    if not exists(file):
        print("error: file not found: " + file)
        exit(4)
    return True

def check_gatk(gatk='~/local/bin/gatk'):
    try:
        gatk = gatk.replace("~", expanduser("~"))
        call([gatk, "--list"],
             stdin=DEVNULL, stdout=DEVNULL, stderr=DEVNULL)
    except FileNotFoundError as e:
        print("check if gatk is installed.")
        exit(134)

    return(gatk)


def check_hisat2_installation():
    out = check_output(["hisat2", "--version"], universal_newlines=True, stderr=subprocess.STDOUT)
    if out.find("hisat2: command not found") > 0:
            print("error: hisat2 not found")
            exit(202)

    # with TemporaryFile('w+t') as e:
    #
    #     call(["hisat2", "--version"], stderr=e)
    #     line = e.readline()
    #
    #     if line.find("hisat2: command not found") >0:
    #         print("error: hisat2 not found")
    #         exit(202)
