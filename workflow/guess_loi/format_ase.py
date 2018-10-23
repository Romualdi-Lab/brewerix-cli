from sys import argv, stdout


def format_ase():
    with open(argv[1], 'r') as fd:
        format_ase_internal(fd, intermediate=stdout)


def format_ase_internal(fd, intermediate, min_coverage=5):
    for idx, line in enumerate(fd):
        tokens = line.rstrip('\n').split('\t')
        if len(tokens) < 7:
            print("less than 7 tokens at line " + str(idx + 1))

        gene_id = '_'.join(tokens[:5])
        ref = tokens[5]
        alt = tokens[6]
        if idx == 0:
            ref_alt = ','.join([ref, alt])
            intermediate.write(gene_id + '\t' + ref_alt + '\n')
        else:
            if int(ref) + int(alt) >= min_coverage:
                ref_alt = ','.join([ref, alt])
                intermediate.write(gene_id + '\t' + ref_alt + '\n')


if __name__ == '__main__':
    format_ase()
