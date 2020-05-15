import unittest
from ncov.parser.qc import is_indel, get_total_variants, is_variant_iupac, is_variant_n, get_qc_data

class TestIndel(unittest.TestCase):
    def test_is_indel_fail(self):
        self.assertFalse(is_indel('A'), 'Should be False, A is not an indel')
    
    def test_is_indel_succeed(self):
        self.assertTrue(is_indel('+A'), 'Should be True, +A is an indel')
    
    def test_count_variants_no_indel(self):
        total_variants = get_total_variants(file='data/sampleA.variants.tsv', indel=False)
        self.assertEqual(total_variants['total_variants'], 9, '9 single nucleotide variant')
        self.assertEqual(total_variants['total_n'], 2, '2 single nucleotide variant')
        self.assertEqual(total_variants['total_iupac'], 4, '4 single nucleotide variant')
    
    def test_count_variants_with_indel(self):
        self.assertEqual(get_total_variants(file='data/sampleA.variants.tsv', indel=True)['total_variants'], 10, '10 variants')
    
    def test_is_variant_n_success(self):
        self.assertTrue(is_variant_n('N'), 'True, variant is N')
        self.assertTrue(is_variant_n('n'), 'True, variant is n')
    
    def test_is_variant_n_fail(self):
        self.assertFalse(is_variant_n('t'), 'False, variant is t')
        self.assertFalse(is_variant_n('G'), 'False, variant is G')

    def test_get_qc_data_pass(self):
        qc_data = get_qc_data(file='data/sampleA.qc.csv')
        self.assertEqual(qc_data['sample_name'], 'sampleA', 'Sample name is: sampleA')
        self.assertEqual(qc_data['pct_covered_bases'], '68.01', 'Percent bases covered: 68.01')
        self.assertEqual(qc_data['qc_pass'], 'FALSE', 'QC pass status: FAIL')


if __name__ == '__main__':
    unittest.main()
