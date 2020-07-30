'''
Suite of tests for test for the ncov.parser.qc module
'''
import unittest
from ncov.parser.qc import *
from ncov.parser.ont import *
from ncov.parser.illumina import *

class TestQc(unittest.TestCase):
    '''
    A unittest class for the QC module
    '''
    def test_is_indel_fail(self):
        '''
        A method to test the is_indel function, should fail this test.
        '''
        self.assertFalse(is_indel('A'), 'Should be False, A is not an indel')

    def test_is_indel_succeed(self):
        '''
        A method to test the is_indel function, should pass this test.
        '''
        self.assertTrue(is_indel('+A'), 'Should be True, +A is an indel')

    def test_count_variants_no_indel(self):
        '''
        A method to test the get_total_variants function without indels.
        '''
        var_file = 'data/illumina/sampleA.variants.tsv'
        fasta_file = 'data/tester.fa'
        total_variants = get_total_variants(file=var_file, reference=fasta_file,
                                            indel=False)
        self.assertEqual(total_variants['total_variants'], 9,
                         '9 single nucleotide variant')

    def test_count_variants_with_indel(self):
        '''
        A method to test the get_total_variants function with indels.
        '''
        var_file = 'data/illumina/sampleA.variants.tsv'
        fasta_file = 'data/tester.fa'
        self.assertEqual(get_total_variants(file=var_file, reference=fasta_file,
                                            indel=True)['total_variants'],
                         10,
                         '10 variants')

    def test_is_variant_n_success(self):
        '''
        A method to test the is_variant_n function, should pass.
        '''
        self.assertTrue(is_variant_n('N'), 'True, variant is N')
        self.assertTrue(is_variant_n('n'), 'True, variant is n')

    def test_is_variant_n_fail(self):
        '''
        A method to test the is_variant_n function, should fail.
        '''
        self.assertFalse(is_variant_n('t'), 'False, variant is t')
        self.assertFalse(is_variant_n('G'), 'False, variant is G')

    def test_get_qc_data_pass(self):
        '''
        A method to test the get_qc_data function
        '''
        qc_data = get_qc_data(file='data/illumina/sampleA.qc.csv')
        self.assertEqual(qc_data['sample_name'], 'sampleA', 'Sample name is: sampleA')
        self.assertEqual(qc_data['pct_covered_bases'], '68.01', 'Percent bases covered: 68.01')
        self.assertEqual(qc_data['qc_pass'], 'FALSE', 'QC pass status: FAIL')

    def test_import_metadata(self):
        '''
        A method to test the import_metadata function.
        '''
        metafile = 'data/metadata.tsv'
        meta_data = import_ct_data(file=metafile)
        self.assertEqual(meta_data['sampleA'], "17.4")
        self.assertEqual(meta_data['sampleB'], "18.4")
        self.assertEqual(len(meta_data.keys()), 2)

    def test_create_qc_summary_line(self):
        '''
        A method to test the create_qc_summary_line function.
        '''
        qcfile = 'data/illumina/sampleA.qc.csv'
        metafile = 'data/metadata.tsv'
        varfile = 'data/illumina/sampleA.variants.tsv'
        covfile = 'data/illumina/sampleA.per_base_coverage.bed'
        fasta = 'data/tester.fa'
        qc_summary = create_qc_summary_line(var_file=varfile,
                                            qc_file=qcfile,
                                            cov_file=covfile,
                                            meta_file=metafile,
                                            indel=True,
                                            fasta=fasta)
        self.assertEqual(qc_summary['sample_name'], 'sampleA', 'Sample name correctly identified \
                         as: sampleA')
        # self.assertEqual(qc_summary['pct_n_bases'], '31.29', 'pct_n_bases is 31.29')
        self.assertEqual(qc_summary['pct_n_bases'], '10.08', 'pct_n_bases is 10.08')
        self.assertEqual(qc_summary['pct_covered_bases'], '68.01', 'pct_covered_bases is 68.01')
        self.assertEqual(qc_summary['total_variants'], 10, 'total_variants is correct')
        self.assertEqual(qc_summary['total_snv'], 9, 'total_snv is correct')
        self.assertEqual(qc_summary['total_indel'], 1, 'total_indel is correct')
        self.assertEqual(qc_summary['total_n'], 13, 'total_n is correct')
        self.assertEqual(qc_summary['total_iupac'], 9, 'total_iupac is correct')
        self.assertEqual(qc_summary['mean_depth'], 679.4, 'mean_depth is correct')
        self.assertEqual(qc_summary['median_depth'], 682, 'median_depth is correct')
        self.assertEqual(qc_summary['ct'], '17.4', 'ct value is correct')
        self.assertEqual(qc_summary['date'], '2020-03-02', 'ct value is correct')
        self.assertEqual(qc_summary['qc_pass'], 'FALSE', 'qc_pass is correct')

    def test_create_qc_summary_line_no_meta(self):
        '''
        A method to test the create_qc_summary_line_no_meta function
        '''
        qcfile = 'data/illumina/sampleA.qc.csv'
        varfile = 'data/illumina/sampleA.variants.tsv'
        covfile = 'data/illumina/sampleA.per_base_coverage.bed'
        fasta = 'data/tester.fa'
        qc_summary = create_qc_summary_line(var_file=varfile,
                                            qc_file=qcfile,
                                            cov_file=covfile,
                                            indel=True,
                                            fasta=fasta)
        self.assertEqual(qc_summary['sample_name'], 'sampleA', \
                                    'Sample name correctly identified as: sampleA')
        # self.assertEqual(qc_summary['pct_n_bases'], '31.29', 'pct_n_bases is 31.29')
        self.assertEqual(qc_summary['pct_n_bases'], '10.08', 'pct_n_bases is 10.08')
        self.assertEqual(qc_summary['pct_covered_bases'], '68.01', \
                                    'pct_covered_bases is 68.01')
        self.assertEqual(qc_summary['total_variants'], 10, 'total_variants is correct')
        self.assertEqual(qc_summary['total_snv'], 9, 'total_snv is correct')
        self.assertEqual(qc_summary['total_indel'], 1, 'total_indel is correct')
        self.assertEqual(qc_summary['total_n'], 13, 'total_n is correct')
        self.assertEqual(qc_summary['total_iupac'], 9, 'total_iupac is correct')
        self.assertEqual(qc_summary['mean_depth'], 679.4, 'mean_depth is correct')
        self.assertEqual(qc_summary['median_depth'], 682, 'median_depth is correct')
        self.assertEqual(qc_summary['ct'], 'NA', 'ct value is correct')
        self.assertEqual(qc_summary['qc_pass'], 'FALSE', 'qc_pass is correct')

    def test_count_iupac_in_fasta(self):
        '''
        A method for testing the count_iupac_in_fasta function.
        '''
        iupac_file = 'data/tester.fa'
        iupac_count = count_iupac_in_fasta(fasta=iupac_file)
        self.assertEqual(iupac_count['total_n'], 13, 'N count is correct')
        self.assertEqual(iupac_count['total_iupac'], 9, 'IUPAC count is correct')

    def test_get_fasta_sequence_length(self):
        '''
        A method to test get_fasta_sequence_length function.
        '''
        fasta = 'data/tester.fa'
        seq_info = get_fasta_sequence_length(fasta=fasta)
        self.assertEqual(seq_info['genome_length'], 129, 'genome_length is correct')

    def test_is_base_masked(self):
        '''
        A method to test the is_base_masked function, test both True and False
        return values.  Note test FASTA file has a genome 129bp of length.
        '''
        fasta = 'data/tester.fa'
        end_pos = get_fasta_sequence_length(fasta=fasta)['genome_length']
        mask_start = 10
        mask_end = 10
        pos_true = 10
        mask_true = is_base_masked(pos=pos_true, end=end_pos,
                                   mask_start=mask_start, mask_end=mask_end)
        self.assertTrue(mask_true, 'position is masked')
        pos_true = 125
        mask_true = is_base_masked(pos=pos_true, end=end_pos,
                                   mask_start=mask_start, mask_end=mask_end)
        self.assertTrue(mask_true, 'position is masked')
        pos_false = 66
        mask_false = is_base_masked(pos=pos_false, end=end_pos,
                                    mask_start=mask_start, mask_end=mask_end)
        self.assertFalse(mask_false, 'position is not masked')

    def test_is_indel_triplet(self):
        '''
        A method to test the is_indel_triplet function.
        '''
        variant_fail = '+AAGGG'
        variant_succeed = '-AAT'
        self.assertFalse(is_indel_triplet(variant_fail), 'variant is not a triplet')
        self.assertTrue(is_indel_triplet(variant_succeed), 'variant is a triplet')

    def test_get_qc_from_consensus_data(self):
        '''
        Validate the output from the get_qc_ont_data
        '''
        fasta = 'data/ont/sampleA.consensus.fasta'
        reference='data/reference.fa'
        consensus_suffix='.consensus.fasta'
        qc_data = get_qc_ont_data(fasta=fasta,
                                       reference=reference,
                                       consensus_suffix=consensus_suffix)
        self.assertEqual(qc_data['sample_name'], 'sampleA', 'valid sample name')
        self.assertEqual(qc_data['pct_n_bases'], '10.08', 'valid pct_n_bases')
        self.assertEqual(qc_data['pct_covered_bases'], 'NA',
                         'valid pct_covered_bases')
        self.assertEqual(qc_data['qc_pass'], 'NA', 'valid qc_pass value')

    def test_get_total_variants_vcf(self):
        '''
        Test the output from the get_total_variants_vcf function.
        '''
        vcf_file = 'data/ont/sampleA.pass.vcf'
        reference = 'data/tester.fa'
        vars = get_total_variants_vcf(file=vcf_file, reference=reference,
                                      mask_start=20, mask_end=20, indel=True)
        self.assertEqual(vars['total_variants'], 9, 'valid number of variants')
        self.assertEqual(vars['total_snv'], 5, 'valid number of snv')
        self.assertEqual(vars['total_indel'], 4, 'valid number of indels')
        self.assertEqual(vars['total_snv_masked'], 3, 'valid number of masked snvs')
        self.assertEqual(vars['total_indel_masked'], 1,
                         'valid number of masked indels')
        self.assertEqual(vars['total_indel_triplet'], 1,
                         'valid number of indel triplets')



if __name__ == '__main__':
    unittest.main()
