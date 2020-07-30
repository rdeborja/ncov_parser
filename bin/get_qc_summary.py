#!/usr/bin/env python
'''
A Python package for summarizing QC data from the ncov-tools pipeline.
'''


import argparse
import sys
import ncov.parser as qc

parser = argparse.ArgumentParser(description="Tool for summarizing QC data")
parser.add_argument('-v', '--variants',
                    help='<sample>.variants.tsv file to process')
parser.add_argument('-e', '--coverage',
                    help='<sample>.per_base_coverage.bed file to process')
parser.add_argument('-c', '--qc', default=None,
                    help='<sample>.qc.csv file to process')
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
parser.add_argument('-n', '--instrument', default='illumina',
                    help='sequencing platform used {illumina, ont}')
if len(sys.argv) == 1:
    parser.print_help(sys.stderr)
    sys.exit(1)
args = parser.parse_args()
qc_line = {}
if args.instrument == 'illumina':
    qc_line.update(qc.get_total_variants(file=args.variants,
                                         reference=args.reference,
                                         indel=args.indel,
                                         mask_start=int(args.mask_start),
                                         mask_end=int(args.mask_end)))
#    qc_line.update(qc.get_qc_data(file=args.qc))
    qc_line.update(qc.get_qc_from_consensus_data(fasta=args.fasta,
                                                 reference=args.reference,
                                                 consensus_suffix='.primertrimmed.consensus.fa'))
elif args.instrument == 'ont':
    qc_line.update(qc.get_total_variants_vcf(file=args.variants,
                                             reference=args.reference,
                                             indel=True,
                                             mask_start=int(args.mask_start),
                                             mask_end=int(args.mask_end)))
    # qc_line.update(qc.get_qc_ont_data(fasta=args.fasta, reference=args.reference,
    #                                   consensus_suffix='.consensus.fa'))
    qc_line.update(qc.get_qc_from_consensus_data(fasta=args.fasta, reference=args.reference,
                                                 consensus_suffix='.consensus.fasta'))
else:
    print("Invalid instrument, exiting...")
    sys.exit(1)

qc_line.update(qc.get_coverage_stats(file=args.coverage))
qc_line.update(qc.count_iupac_in_fasta(fasta=args.fasta))
try:
    meta_data = qc.import_metadata(file=args.meta)
    qc_line['ct'] = meta_data[qc_line['sample_name']]['ct']
    qc_line['date'] = meta_data[qc_line['sample_name']]['date']
except:
    qc_line['ct'] = 'NA'
    qc_line['date'] = 'NA'

qc.write_qc_summary_header()
qc.write_qc_summary(summary=qc_line)
