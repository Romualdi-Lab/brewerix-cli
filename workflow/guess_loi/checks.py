from os.path import exists
from shutil import which


def check_file_exists(file):
    if not exists(file):
        print("error: file not found: " + file)
        exit(4)
    return True


def check_command_availability(commands):
    for cmd in commands:
        if which(cmd) is None:
            exit("program not installed: %s" % cmd)


def vcf_index_exits(file):
    return exists(file)