'''
A set of functions used in the Illumina pipeline.
'''

import re
import statistics
import glob
import csv
import vcf
import os
from Bio import SeqIO
from .qc import *

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


def get_total_variants(file, reference, mask_start=100, mask_end=50,
                       indel=False):
    '''
    A function that parses the iVar variants file and returns the total number
    of variants.

    Arguments:
        * file:         a string containing the filename and path to the
                        <sample>.variants.tsv file
        * reference:    full path to the reference FASTA file
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
    total_variants = 0
    try:
        for record in SeqIO.parse(reference, 'fasta'):
            genome_length = len(str(record.seq))
    except:
        genome_length = 0
    with open(file) as file_p:
        variant_reader = csv.DictReader(file_p, delimiter='\t')
        for data in variant_reader:
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
        total_variants = counter_snv + counter_indel
    file_p.close()
    return {'total_variants' : total_variants,
            'total_snv' : counter_snv,
            'total_indel' : counter_indel,
            'total_snv_masked' : counter_snv_masked,
            'total_indel_masked' : counter_indel_masked,
            'total_indel_triplet' : counter_indel_triplet,
            'genome_length' : genome_length}


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
            * mean_depth
            * median_depth
            * ct
            * date
            * qc_pass
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

    # import the ct and collection date from the metadata file
    try:
        #meta_data = import_ct_data(file=meta_file)
        meta_data = import_metadata(file=meta_file)
        summary['ct'] = meta_data[summary['sample_name']]['ct']
        summary['date'] = meta_data[summary['sample_name']]['date']
    except:
        summary['ct'] = 'NA'
        summary['date'] = 'NA'
    return summary
