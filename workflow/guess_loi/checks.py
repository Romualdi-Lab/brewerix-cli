from subprocess import call

from os.path import exists


def check_file_exists(file):
    if not exists(file):
        print("error: file not found: " + file)

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
    with TemporaryFile() as e:
        call(["hisat2", "--version"], stderr=e)

    line = e.readline()
    if line.find("hisat2: command not found") >0:
        print("error: hisat2 not found")
        exit(202)
