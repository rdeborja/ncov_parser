#!/usr/bin/env python
'''
A Python package for summarizing QC data from the ncov-tools pipeline.
'''


import argparse
import sys
from ncov.parser.qc import create_qc_summary_line, write_qc_summary, \
write_qc_summary_header, import_ct_data


parser = argparse.ArgumentParser(description="Tool for summarizing QC data")
parser.add_argument('-c', '--qc', help='<sample>.qc.csv file to process')
parser.add_argument('-v', '--variants', \
        help='<sample>.variants.tsv file to process')
parser.add_argument('-e', '--coverage', \
        help='<sample>.per_base_coverage.bed file to process')
parser.add_argument('-i', '--indel', action='store_true', \
        help='flag to determine whether to count indels')
parser.add_argument('-m', '--meta', \
        help='full path to the metadata YAML file')

if len(sys.argv) == 1:
    parser.print_help(sys.stderr)
    sys.exit(1)

args = parser.parse_args()

ct_data = import_ct_data(file=args.meta)
qc_line = create_qc_summary_line(
        var_file=args.variants,
        qc_file=args.qc,
        cov_file=args.coverage,
        meta_file=args.meta,
        indel=args.indel)

write_qc_summary_header()
write_qc_summary(summary=qc_line)
