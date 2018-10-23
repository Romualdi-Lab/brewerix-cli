from re import sub, search


def guess_sample_name(fq, paired=False):
    if paired:
        return sub(r'_1\.fq(\.gz)*$', '', fq)
    return sub(r'\.fq(\.gz)*$', '', fq)


def check_paired_end_nomenclature(fq):
    if search(r'_R[12]', fq):
        print("error: paired end nomenclature format: use _1/_2 not _R1/_R2")
        exit(202)
