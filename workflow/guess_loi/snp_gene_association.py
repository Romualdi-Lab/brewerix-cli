import re
from sys import argv


def snp2gene():
    vcf_file = argv[1]
    snp2gene_association_from_vcf(vcf_file)

def snp2gene_association_from_vcf(vcf_file):
    with open(vcf_file, 'r') as fd:
        for idx, line in enumerate(fd):
            if line[0] == '#':
                continue

            tokens = line.rstrip('\n').split('\t')
            if tokens[7].find('GENEINFO=') == -1:
                print("No gene information fuond at line " + str(idx+1))
                exit(559)
            base_id = '_'.join(tokens[:6])
            p = re.compile(r'GENEINFO=([^;]+);')
            gene_infos = p.findall(tokens[7])

            for gene_info in gene_infos[0].split('|'):
                entrez, symbol = gene_info.split(':')
                print(base_id + '\t' + symbol + '\t' + entrez)

if __name__ == '__main__':
    snp2gene()
