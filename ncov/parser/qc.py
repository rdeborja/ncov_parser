'''
A parser for the qc metrics file generated by the COG-UK Nextflow pipeline.  The
generated file contains a header line and a data line with pre-defined columns.
'''

import re
import statistics
import glob
import csv
from Bio import SeqIO


def get_qc_data(file):
    '''
    A function to parse the COG-UK QC file and returns a data structure with the
    QC results.

    Arguments:
        * file (str): a string containing the path and filename of the <sample>.qc.csv file

    Return Value:
        * dict: returns a dictionary with keys "sample_name", "pct_covered_bases", "qc_pass"
    '''
    with open(file) as file_p:
        qc_reader = csv.DictReader(file_p, delimiter=',')
        for line in qc_reader:
            sample_name = line['sample_name']
            pct_n_bases = line['pct_N_bases']
            pct_covered_bases = line['pct_covered_bases']
            qc_pass = line['qc_pass']
    file_p.close()
    return {'sample_name' : sample_name,
            'pct_n_bases' : pct_n_bases,
            'pct_covered_bases' : pct_covered_bases,
            'qc_pass' : qc_pass}


def get_total_variants(file, reference, start=1, mask_start=100, mask_end=50,
                       indel=False):
    '''
    A function that parses the iVar variants file and returns the total number
    of variants.

    Arguments:
        * file:         a string containing the filename and path to the
                        <sample>.variants.tsv file
        * reference:    full path to the reference FASTA file
        * start:        genome start position (default: 1)
        * mask_start:   bases to mask at beginning of genome (default: 100)
        * mask_end:     bases to mask at end of genome (default: 50)
        * indel:        a boolean to determine whether to process indels

    Returns:
        Function returns a dictionary containing the following keys:
            * total_variants:       total number of variants in the file
            * total_snv:            total number of SNV variants in the file
            * total_indel:          total number of indel variants in the file
            * total_snv_masked:     total number of masked SNVs
            * total_indel_masked:   total number of masked indels
            * genome_length:        length of genome sequence in FASTA file
    '''
    counter = 0
    counter_snv = 0
    counter_indel = 0
    counter_snv_masked = 0
    counter_indel_masked = 0
    genome_length = 0
    counter_indel_triplet = 0
    try:
        for record in SeqIO.parse(reference, 'fasta'):
            genome_length = len(str(record.seq))
    except:
        genome_length = 0
    with open(file) as file_p:
        variant_reader = csv.DictReader(file_p, delimiter='\t')
        for data in variant_reader:
        # for line in file_p:
            # if re.match("^REGION\tPOS\tREF", line):
                # skip to the next line if header encountered
                # continue
            # check if the variant is an indel and the option for counting
            # indels
            # data = line.split("\t")
            base_masked = is_base_masked(pos=int(data['POS']),
                                         end=genome_length,
                                         mask_start=mask_start,
                                         mask_end=mask_end)
            if indel:
                if len(str(data['ALT'])) > 1:
                    counter += 1
                    counter_indel += 1
                    if is_indel_triplet(data['ALT']):
                        counter_indel_triplet += 1
                    if genome_length > 0:
                        if base_masked:
                            counter_indel_masked += 1
                elif len(str(data['ALT'])) == 1:
                    counter += 1
                    counter_snv += 1
                    if genome_length > 0:
                        if base_masked:
                            counter_snv_masked += 1
            elif not indel:
                if len(str(data['ALT'])) == 1:
                    counter += 1
                    counter_snv += 1
                    if genome_length > 0:
                        if base_masked:
                            counter_snv_masked += 1
                else:
                    continue
    file_p.close()
    return {'total_variants' : counter,
            'total_snv' : counter_snv,
            'total_indel' : counter_indel,
            'total_snv_masked' : counter_snv_masked,
            'total_indel_masked' : counter_indel_masked,
            'total_indel_triplet' : counter_indel_triplet,
            'genome_length' : genome_length}


def get_variant_dictionary():
    '''
    Create a dictionary of variants per sample.
    '''



def import_ct_data(file, sample_id='sample', ct_id='ct'):
    '''
    Obtain the name of the metadata YAML file and import the data.  This
    assumes a header exists in the file.

    Arguments:
        * file:         full path to the metadata YAML file
        * sample_id:    the column label representing the sample name
                        (default: 'sample')
        * ct_id:        the column label representing the ct value
                        (default: 'ct')

    Return Value:
        The function returns a dictionary with {"sample" : "ct"}
    '''
    try:
        with open(file) as file_p:
            data = {}
            ct_reader = csv.DictReader(file_p, delimiter='\t')
            for line in ct_reader:
                data[line[sample_id]] = line[ct_id]
        return data
    except:
        return {'ct' : 'NA'}


def is_variant_n(variant):
    '''
    A function to determine whether the mutation is N

    Arguments:
        * variant: a string representing the variant

    Return Value:
        Function returns a boolean
    '''
    variant = str(variant).upper()
    return re.search('[Nn]', variant)


def is_variant_iupac(variant):
    '''
    A function to determine whether a variant is an IUPAC code, note that we
    are treating N as a distinct value.

    Arguments:
        * variant: a string reprenting the variant

    Return Value:
        Function returns a boolean
    '''
    variant = str(variant).upper()
    iupac_codes = '[RYSWKMBDHV]'
    return re.search(iupac_codes, variant)


def count_iupac_in_fasta(fasta):
    '''
    Count the number of IUPAC occurrences, not including [Nn], in the consensus
    FASTA file.

    Arguments:
        * fasta:    FASTA reference file

    Return Value:
        The function returns an integer representing the number of IUPAC code
        occurrences in the FASTA file.  Note that Ns were not considered.
    '''
    iupac_codes = '[RYSWKMBDHV]'
    try:
        for record in SeqIO.parse(fasta, 'fasta'):
            iupac_count = re.subn(iupac_codes, repl='', string=str(record.seq))[1]
            n_count = re.subn('[Nn]', repl='', string=str(record.seq))[1]
            consensus_length = len(str(record.seq))
        return {'total_n' : n_count,
                'total_iupac' : iupac_count,
                'consensus_length' : consensus_length}
    except:
        return {'total_n' : 'NA',
                'total_iupac' : 'NA',
                'consensus_length' : 'NA'}


def get_fasta_sequence_length(fasta):
    '''
    Return the length of the genome in the FASTA sequence file.

    Arguments:
        * fasta: full path to the FASTA reference genome file

    Return Value:
        Returns an integer representing the lenght of the reference genome
        from the FASTA file.
    '''
    try:
        for record in SeqIO.parse(fasta, 'fasta'):
            seq_length = len(record.seq)
        return {'genome_length' : seq_length}
    except:
        return {'genome_length' : 'NA'}


def is_base_masked(pos, end, start=1, mask_start=100, mask_end=50):
    '''
    Check if the base is masked within the beginning and end of the sequence.
    The ends are considered difficult to sequence and will be masked from
    downstream analysis.

    Arguments:
        * pos:          position of the base
        * start:        start position of the genome (default: 1)
        * end:          end position of the genome
        * mask_start:   bases to mask from beginning of sequence (default: 100)
        * mask_end:     bases to mask at end of sequence (default: 50)

    Return Value:
        Function returns a boolean
    '''
    return int(pos) < (int(start) + int(mask_start)) or int(pos) > (int(end) - int(mask_end))


def is_indel(variant):
    '''
    Check whether the variant from the <sample>.variants.tsv file is an indel.
    Note that indels will have a +/- in the ALT column of the file.

    Arguments:
        * variant: a string representing the variant

    Return Value:
        Function returns a boolean
    '''
    return len(variant) > 1


def is_indel_triplet(variant, size=3):
    '''
    Check whether the indel is 3bp in size.  We will be using this to infer potential
    frameshift indels.

    Arguments:
        * variants: a string represent the variant
        * size:     size of a codon indicating a frameshift

    Return Value:
        Function returns a boolean
    '''
    # remove the leading +/- from the variant
    variant = re.sub('^[+-]', '', variant)
    if len(variant) > 0:
        return not len(variant) % size


def get_coverage_stats(file):
    '''
    A function to calculate the depth of coverage across the genome from the
    bedtools <sample>.per_base_coverage.bed file.

    Arguments:
        * file: a string containing the filename and path to the
                <sample>.per_sample_coverage.bed file

    Return Value:
        Function returns a dictionary with the following keys:
            * mean_depth
            * median_depth
    '''
    depth = []
    with open(file) as file_p:
        for line in file_p:
            if re.match("^reference_name\tstart\tend", line):
                # skip to the next line if header encountered
                continue
            line = line.strip()
            data = line.split("\t")
            depth.append(int(data[7]))
    file_p.close()
    mean_depth = round(statistics.mean(depth), 1)
    median_depth = round(statistics.median(depth), 1)
    return {"mean_depth" : mean_depth, "median_depth" : median_depth}


def create_qc_summary_line(var_file, qc_file, cov_file, meta_file=None,
                           indel=True, fasta=None, mask_start=100,
                           mask_end=50, reference=None):
    '''
    A function that aggregates the different QC data into a single sample
    dictionary entry.

    Arguments:
        * var_file:     full path to the <sample>.variants.tsv file
        * qc_file:      full path to the <sample>.qc.csv file
        * cov_file:     full path to the <sample>.per_base_coverage.bed file
        * meta_file:    full path to the 'metadata.tsv' file
        * fasta:        full path to the <sample>.consensus.fa file
        * reference:    full path to the reference FASTA genome file
        * indel:        boolean to determine whether to use indels in variant
                        count (default: True)

    Return Value:
        Return an aggregate dictionary containing the following keys:
            * total_variants
            * total_snv
            * total_indel
            * total_n
            * total_iupac
            * consensus_length
            * sample_name
            * pct_n_bases
            * pct_covered_bases
            * qc_pass
            * mean_depth
            * median_depth
            * ct
    '''
    summary = {}
    summary.update(get_total_variants(file=var_file,
                                      indel=indel,
                                      reference=reference,
                                      mask_start=mask_start,
                                      mask_end=mask_end))
    summary.update(get_qc_data(file=qc_file))
    summary.update(get_coverage_stats(file=cov_file))
    summary.update(count_iupac_in_fasta(fasta=fasta))

    # import the ct data from the metadata file and add the ct value to the summary dictionary
    try:
        meta_data = import_ct_data(file=meta_file)
        summary['ct'] = meta_data[summary['sample_name']]
    except:
        summary['ct'] = 'NA'
    return summary


def write_qc_summary(summary):
    '''
    A function to write the QC data line to output in the order:
    * sample name
    # % n bases
    * % bases covered
    * total mutations
    * total snv
    * total snv masked
    * total indel
    * total indel masked
    * total N mutations
    * total IUPAC mutations
    * mean sequence depth
    * median sequence depth
    * iVar QC pass

    Arguments:
        * summary:  dictionary containing the keys sample_name, pct_n_bases,
                    pct_covered_bases, total_variants, total_snv, total_indel,
                    total_n, total_iupac, mean_depth, median_depth, ct, qc_pass

    Return Value:
        None
    '''
    summary_line = '\t'.join([
        summary['sample_name'],
        str(summary['pct_n_bases']),
        str(summary['pct_covered_bases']),
        str(summary['total_variants']),
        str(summary['total_snv']),
        str(summary['total_snv_masked']),
        str(summary['total_indel']),
        str(summary['total_indel_masked']),
        str(summary['total_indel_triplet']),
        str(summary['total_n']),
        str(summary['total_iupac']),
        str(summary['mean_depth']),
        str(summary['median_depth']),
        str(summary['ct']),
        str(summary['qc_pass'])])
    print(summary_line)


def write_qc_summary_header(header=['sample_name',
                                    'pct_n_bases',
                                    'pct_covered_bases',
                                    'total_variants',
                                    'total_snv',
                                    'total_snv_masked',
                                    'total_indel',
                                    'total_indel_masked',
                                    'total_indel_triplet',
                                    'total_n',
                                    'total_iupac',
                                    'mean_depth',
                                    'median_depth',
                                    'ct',
                                    'qc_pass']):
    '''
    Write the header for the QC summary data

    Arguments:
        * header: a list containing the column header

    Return Value:
        None
    '''
    print('\t'.join(header))


def collect_qc_summary_data(path, pattern='.summary.qc.tsv'):
    '''
    An aggregation function to collect individual sample based QC summary data
    and create a single file with all samples.

    Arguments:
        * path:     full path to the <sample>.summary.qc.tsv files
        * pattern:  file pattern for the sample files (default: .summary.qc.tsv)

    Return Value:
        data: a list containing the summary line data
    '''
    files = glob.glob(path + "/*" + pattern)
    data = []
    for file in files:
        with open(file) as file_p:
            for line in file_p:
                # skip the header
                if re.match("^sample_name\tpct_n_bases\tpct_covered_bases", line):
                    continue
                data.append(line.rstrip())
    return data
