'''

'''

#!/usr/bin/env python

import argparse
import sys
from ncov.parser.qc import collect_qc_summary_data, write_qc_summary_header

parser = argparse.ArgumentParser(description="Tool for aggregating sample QC data")
parser.add_argument('-p', '--path', \
    help='directory to search for <sample>.summary.qc.tsv files')

if len(sys.argv) == 1:
    parser.print_help(sys.stderr)
    sys.exit(1)

args = parser.parse_args()

summary_data = collect_qc_summary_data(dir=args.path)

write_qc_summary_header()
for summary_line in summary_data:
    print(summary_line)
