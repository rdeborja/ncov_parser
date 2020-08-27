# ncov_parser

# THIS PACKAGE HAS BEEN DEPRECATED AND WILL NO LONGER BE UPDATED
PLEASE SEE https://github.com/jts/ncov-tools for an updated version of the package.


[![Build Status](https://travis-ci.com/rdeborja/ncov_parser.svg?branch=master)](https://travis-ci.com/rdeborja/ncov_parser) [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

The `ncov_parser` package provides a suite of tools to parse the files generated
in the Nextflow workflow and provide a QC summary file.  The package requires
several files including:
* <sample>.variants.tsv
* <sample>.qc.csv
* <sample>.per_base_coverage.bed
* <sample>.primertrimmed.consensus.fa
* <reference genome>.fa

An optional metadata file with `ct` values can be included.

In addition, `bedtools` should be run to generate a
`<sample>.per_base_coverage.bed` file to generate mean and median depth of
coverage statistics.


## Installation
After downloading the repository, the package can be installed using `pip`:
```
git clone git@github.com:rdeborja/ncov_parser.git
cd ncov_parser
pip install .
```


## Usage
The library consists of several functions that can be imported.
```
from ncov.parser.qc import *
```
Similarly, you can import only those functions of interesting, this can include:
```
get_qc_data
get_total_variants
import_ct_data
is_variant_n
is_variant_iupac
count_iupac_in_fasta
get_fasta_sequence_length
is_base_masked
is_indel
is_indel_triplet
get_coverage_stats
create_qc_summary_line
write_qc_summary
write_qc_summary_header
collect_qc_summary_data
```

### Top levels scripts
In the `bin` directory, several wrapper scripts exist to assist in generating
QC metrics.

To create sample level summary qc files, use the `get_qc_summary.py` script:
```
get_qc_summary.py --qc <sample>.qc.csv --variants <sample>.variants.tsv
--coverage <sample>.per_base_coverage.bed --meta <metadata>.tsv
--fasta <sample>.primertrimmed.consensus.fa --reference <reference genome>.fa
[--indel] [--mask_start 100] [--mask_end 50]
```

Note the `--indel` flag should only be present if indels will be used in the
calculation of variants.

Once this is complete, we can use the `collect_qc_summary.py` script to
aggregate the sample level summary files into a single run tab-separate file.
```
collect_qc_summary.py --path <path to sample.summary.qc.tsv files>
```

Note that this tool has been used in conjunction with the [@jts `ncov-tools`](https://github.com/jts/ncov-tools)
suite of tools. 

## License
MIT
