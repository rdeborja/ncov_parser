# ncov_parser

The `ncov_parser` package provides a suite of tools to parse the files generated
in the Nextflow workflow and provide a QC summary file.  The package requires
several files including:
* <sample>.variants.tsv
* <sample>.qc.csv

In addition, `bedtools` should be run to generate a
`<sample>.per_base_coverage.bed` file to generate mean and median depth of
coverage statistics.


# Installation
After downloading the repository, the package can be installed using `pip`:
```
git clone git@github.com:rdeborja/ncov_parser.git
cd ncov_parser
pip install .
```


# Usage
The library consists of several functions that can be imported.
```
from ncov.parser.qc import get_qc_data, get_total_variants, get_coverage_stats,
create_qc_summary_line, write_qc_summary, collect_qc_summary_data
```


# License
MIT

