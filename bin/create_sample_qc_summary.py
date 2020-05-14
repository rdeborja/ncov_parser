#!/usr/bin/env python
'''
A Python package for summarizing QC data from the ncov-tools pipeline.
'''


import argparse
import sys
from ncov.parser.qc import create_qc_summary_line, write_qc_summary, \
write_qc_summary_header


parser = argparse.ArgumentParser(description="Tool for summarizing QC data")
parser.add_argument('-s', '--sample', help='sample name used as a prefix to all files')
parser.add_argument('-c', '--qc_dir', help='<sample>.qc.csv file to process')
parser.add_argument('-v', '--variants_dir', \
        help='<sample>.variants.tsv file to process')
parser.add_argument('-e', '--coverage_dir', \
        help='<sample>.per_base_coverage.bed file to process')
parser.add_argument('-i', '--indel', action='store_true', \
        help='flag to determine whether to count indels')

if len(sys.argv) == 1:
    parser.print_help(sys.stderr)
    sys.exit(1)

args = parser.parse_args()

qc_line = create_qc_summary_line(
        var_file=args.variants_dir + '/' + args.sample + '.variants.tsv',
        qc_file=args.qc_dir + '/' + args.sample + '.qc.csv',
        cov_file=args.coverage_dir + '/' + args.sample + '.per_base_coverage.bed',
        indel=args.indel)

write_qc_summary_header()
write_qc_summary(summary=qc_line)
