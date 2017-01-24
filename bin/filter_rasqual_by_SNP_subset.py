#!/usr/bin/env python3
# Python/3.5.1-foss-2015b
import argparse
import sys
import logging

def main(args):
    SNPs = read_snps_from_bed(args.SNP_bed)
    logging.info('Loaded {} SNPs in filter.'.format(len(SNPs)))
    logging.info('SNP example: "{}"'.format(list(SNPs)[0]))
    for line in filter_SNPs(args.input, SNPs):
        args.output.write(line)

def filter_SNPs(f, SNPs):
    for line in f:
        try:
            snp = line.split('\t')[1]
        except IndexError:
            logging.warning('Unexpected input: {}'.format(line))
            continue
        if snp in SNPs:
            yield line

def read_snps_from_bed(f):
    return {line.split('\t')[-1].strip() for line in f}

if __name__ == '__main__':
    p = argparse.ArgumentParser(description='Filter rasqual output by a set of SNPs of interest')
    p.add_argument('SNP_bed', type=argparse.FileType('r'), 
           help='Bed file with SNP names in 4th column. Only these SNPs are kept from rasqual output.')
    p.add_argument('-i', '--input', default=sys.stdin, type=argparse.FileType('r'),
           help = 'defaults to stdin')
    p.add_argument('-o', '--output', default=sys.stdout, type=argparse.FileType('w'),
    help = 'defaults to stdout')
    args = p.parse_args()
    logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.INFO)
    main(args)
