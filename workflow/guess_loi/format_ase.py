from sys import argv, stdout


def format_ase():
    with open(argv[1], 'r') as fd:
        format_ase_internal(fd, intermediate=stdout)


def format_ase_internal(fd, intermediate):
    for idx, line in enumerate(fd):
        tokens = line.rstrip('\n').split('\t')
        if len(tokens) < 7:
            print("less than 7 tokens at line " + str(idx + 1))

        id = '_'.join(tokens[:5])
        ref = tokens[5]
        alt = tokens[6]
        ref_alt = ','.join([ref, alt])
        intermediate.write(id + '\t' + ref_alt + '\n')


if __name__ == '__main__':
    format_ase()
