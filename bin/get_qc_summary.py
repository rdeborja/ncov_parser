#!/usr/bin/env python
'''
A Python package for summarizing QC data from the ncov-tools pipeline.
'''


import argparse
import sys
from ncov.parser.qc import create_qc_summary_line, write_qc_summary, \
        write_qc_summary_header, import_ct_data, get_qc_data, \
        get_coverage_stats, get_total_variants, count_iupac_in_fasta


parser = argparse.ArgumentParser(description="Tool for summarizing QC data")
parser.add_argument('-c', '--qc', help='<sample>.qc.csv file to process')
parser.add_argument('-v', '--variants',
                    help='<sample>.variants.tsv file to process')
parser.add_argument('-e', '--coverage',
                    help='<sample>.per_base_coverage.bed file to process')
parser.add_argument('-i', '--indel', action='store_true',
                    help='flag to determine whether to count indels')
parser.add_argument('-m', '--meta', default=None,
                    help='full path to the metadata file')
parser.add_argument('-f', '--fasta', default=None,
                    help='full path to the FASTA consensus file')
parser.add_argument('-r', '--reference', default=None,
                    help='full path to the reference FASTA file')
parser.add_argument('--mask_start', default=100,
                    help='number of bases to mask at start of genome')
parser.add_argument('--mask_end', default=50,
                    help='number of bases to mask at end of genome')
if len(sys.argv) == 1:
    parser.print_help(sys.stderr)
    sys.exit(1)
args = parser.parse_args()
ct_data = import_ct_data(file=args.meta)
qc_line = {}
qc_line.update(get_total_variants(file=args.variants,
                                  reference=args.reference,
                                  indel=args.indel,
                                  start=1,
                                  mask_start=int(args.mask_start),
                                  mask_end=int(args.mask_end)))
qc_line.update(get_qc_data(file=args.qc))
qc_line.update(get_coverage_stats(file=args.coverage))
qc_line.update(count_iupac_in_fasta(fasta=args.fasta))
try:
    meta_data = import_ct_data(file=args.meta)
    qc_line['ct'] = meta_data[qc_line['sample_name']]
except:
    qc_line['ct'] = 'NA'

write_qc_summary_header()
write_qc_summary(summary=qc_line)
