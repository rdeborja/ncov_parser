'''

'''

import os
import sys
import csv
from .qc import *


def get_qc_ont_data(fasta, reference, consensus_suffix='.consensus.fa'):
    '''
    In the current pipeline, the <sample>.qc.csv file is not generated for
    nanopore runs.  We can generate the pct_n_bases and the sample_name and
    leave the remaining fields as NA
    '''
    genome_length = 0
    qc_data = {}
    pct_n_bases = ''
    sample_name = re.sub(consensus_suffix, '', os.path.basename(fasta))

    try:
        for record in SeqIO.parse(reference, 'fasta'):
            genome_length = len(str(record.seq))
        qc_data.update(count_iupac_in_fasta(fasta=fasta))
        pct_n_bases = float(qc_data['total_n']) / genome_length * 100
        pct_n_bases = str(round(pct_n_bases, 2))
    except:
        pct_n_bases = 'NA'
    return {'sample_name': sample_name,
            'pct_n_bases' : pct_n_bases,
            'pct_covered_bases' : 'NA',
            'qc_pass' : 'NA'}

def create_vcf_variant_id(var):
    '''
    Create a variant ID consisting of position, reference allele,
    alternate allele as a concatenated string.  Note that the alternate
    allele is provided as a list and needs to be set as a string
    and merged.

    Arguments:
        * var: a VCF record object from pyvcf

    Return Value:
        Returns a string containing a concatenation of variant
        object identifiers.
    '''
    var_alt = []
    for _var in var.ALT:
        var_alt.append(str(_var))
    var_alt_str = ''.join(var_alt)
    return '_'.join([str(var.POS), var.REF, var_alt_str])


def get_total_variants_vcf(file, reference, start=1, mask_start=100,
                           mask_end=50, indel=True):
    '''
    Get the total variants, total snv, total indels, total masked SNVs,
    total masked indels and the genome length from the VCF file.
    '''
    counter = 0
    counter_snv = 0
    counter_indel = 0
    counter_snv_masked = 0
    counter_indel_masked = 0
    genome_length = 0
    counter_indel_triplet = 0
    var_dict = {}
    try:
        for record in SeqIO.parse(reference, 'fasta'):
            genome_length = len(str(record.seq))
    except:
        genome_length = 0

    vcf_reader = vcf.Reader(filename=file)
    for var in vcf_reader:
        # create a unique variant identifier for de-duplication
        var_id = create_vcf_variant_id(var=var)
        if var_id in var_dict:
            continue
        else:
            var_dict[var_id] = 1

        # check if variant position is masked
        base_masked = is_base_masked(pos=int(var.POS),
                                     end=genome_length,
                                     mask_start=mask_start,
                                     mask_end=mask_end)
        # keep a count of total variants
        if indel:
            counter += 1
            if var.is_indel:
                counter_indel += 1
                if is_indel_triplet(re.sub('^.', '', str(var.ALT[0]))):
                    counter_indel_triplet += 1
                elif is_indel_triplet(re.sub('^.', '', str(var.REF))):
                    counter_indel_triplet += 1
                if genome_length > 0:
                    if base_masked:
                        counter_indel_masked += 1
            if var.is_snp:
                counter_snv += 1
                if base_masked:
                    counter_snv_masked += 1
        elif not indel:
            if var.is_snp:
                counter += 1
                counter_snv += 1
                if base_masked:
                    counter_snv_masked += 1
            if var.is_indel:
                continue
        else:
            continue
    return {'total_variants' : counter,
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
