'''
A Python module for handling the metadata.tsv file for the nCoV project.
'''

import csv
import datetime

class Metadata:
    '''
    A class to model the metadata.tsv file.
    '''

    def __init__(self, file, delimiter='\t', sample_id='sample', ct_id='ct',
                 date_id='date'):
        '''
        Class initializer, requires the metadata.tsv file.
        '''
        self.file = file
        self.delimiter = delimiter
        try:
            with open(self.file) as file_p:
                data = {}
                print("Getting data...")
                meta_reader = csv.DictReader(file_p, delimiter=self.delimiter)
                for line in meta_reader:
                    data[line[sample_id]] = {'ct': line[ct_id],
                                             'date' : line[date_id]}
                self.metadata = data
        except:
            self.metadata = None


    def get_sample_metadata(self, sample):
        '''
        Lookup the metadata after importing and return the ct and collection
        date. 
        '''
        if sample in self.metadata:
            return self.metadata[sample]
        else:
            return {'ct' : 'NA', date : 'NA'}


    def get_collection_date(self, date_str, date_fmt='%d%b%Y', output_fmt='%Y-%m-%d'):
        '''
        Format the collection date to the standard format <YEAR>-<MO>-<DA>
        '''
        return datetime.datetime.strptime(date_str, date_fmt).strftime(output_fmt)
