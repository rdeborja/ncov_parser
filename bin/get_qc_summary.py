#!/usr/bin/env python
'''
A Python package for summarizing QC data from the ncov-tools pipeline.
'''


import argparse
import sys
from ncov.parser.qc import create_qc_summary_line, write_qc_summary, \
write_qc_summary_header


parser = argparse.ArgumentParser(description="Tool for summarizing QC data")
parser.add_argument('-c', '--qc', help='<sample>.qc.csv file to process')
parser.add_argument('-v', '--variants', \
        help='<sample>.variants.tsv file to process')
parser.add_argument('-d', '--depth', \
        help='<sample>.per_base_coverage.bed file to process')

if len(sys.argv) == 1:
    parser.print_help(sys.stderr)
    sys.exit(1)

args = parser.parse_args()

qc_line = create_qc_summary_line(
        var_file=args.variants,
        qc_file=args.qc,
        cov_file=args.depth)

write_qc_summary_header()
write_qc_summary(summary=qc_line)
