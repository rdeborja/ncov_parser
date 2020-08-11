'''

'''

import re
import statistics
import glob
import csv
from Bio import SeqIO

class Variants:
    '''
    The Variants class to handle the variants.tsv file.
    '''

    def __init__(self, file):
        '''
        Initialize the object with the variants.tsv file
        '''
        self.file = file


    def get_total_variants(self, indel=True):
        '''
        A method that parses the iVar variants file and returns the total
        number of variants.

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
                * total_snv:            total number of SNV variants in the
                                        file
                * total_indel:          total number of indel variants in the
                                        file
                * total_snv_masked:     total number of masked SNVs
                * total_indel_masked:   total number of masked indels
                * genome_length:        length of genome sequence in FASTA file
        '''
        counter = 0
        counter_snv = 0
        counter_indel = 0
        counter_indel_triplet = 0

        with open(self.file) as file_p:
            variant_reader = csv.DictReader(file_p, delimiter='\t')
            for data in variant_reader:
                if indel:
                    if len(str(data['ALT'])) > 1:
                        counter += 1
                        counter_indel += 1
                        if self.is_indel_triplet(data['ALT']):
                            counter_indel_triplet += 1
                    elif len(str(data['ALT'])) == 1:
                        counter += 1
                        counter_snv += 1
                elif not indel:
                    if len(str(data['ALT'])) == 1:
                        counter += 1
                        counter_snv += 1
                    else:
                        continue
        file_p.close()
        return {'num_variants_snvs' : counter_snv,
                'num_variants_indel' : counter_indel,
                'num_variants_indel_triplet' : counter_indel_triplet}


    def is_indel_triplet(self, variant, size=3):
        '''
        Check whether the indel is a multiple of 3bp in size.  We will be using
        this to infer potential
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
