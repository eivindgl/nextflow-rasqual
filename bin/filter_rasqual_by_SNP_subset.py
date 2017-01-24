#!/usr/bin/env python3
# Python/3.5.1-foss-2015b
import argparse
import sys

def main(args):
   print(args) 


if __name__ == '__main__':
   p = argparse.ArgumentParser(description='Filter rasqual output by a set of SNPs of interest')
   p.add_argument('SNP_bed', type=argparse.FileType('r'), 
           help='Bed file with SNP names in 4th column. Only these SNPs are kept from rasqual output.')
   p.add_argument('-i', '--input', default=sys.stdin, type=argparse.FileType('r'),
           help = 'defaults to stdin')
   p.add_argument('-o', '--output', default=sys.stdout, type=argparse.FileType('w'),
    help = 'defaults to stdout')
   args = p.parse_args()
   main(args)
