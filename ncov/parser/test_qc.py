import unittest
from ncov.parser.qc import is_indel, get_total_variants, is_variant_iupac, \
    is_variant_n, get_qc_data, import_ct_data, create_qc_summary_line, \
    count_iupac_in_fasta

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
        self.assertEqual(get_total_variants(file='data/sampleA.variants.tsv',indel=True)['total_variants'], 10, '10 variants')

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

    def test_import_metadata(self):
        metafile = 'data/metadata.tsv'
        meta_data = import_ct_data(file=metafile)
        self.assertEqual(meta_data['sampleA'], "17.4")
        self.assertEqual(meta_data['sampleB'], "18.4")
        self.assertEqual(len(meta_data.keys()), 2)

    def test_create_qc_summary_line(self):
        qcfile = 'data/sampleA.qc.csv'
        metafile = 'data/metadata.tsv'
        varfile = 'data/sampleA.variants.tsv'
        covfile = 'data/sampleA.per_base_coverage.bed'
        qc_summary = create_qc_summary_line(var_file=varfile,
                                            qc_file=qcfile,
                                            cov_file=covfile,
                                            meta_file=metafile,
                                            indel=True)
        self.assertEqual(qc_summary['sample_name'], 'sampleA', 'Sample name correctly identified as: sampleA')
        self.assertEqual(qc_summary['pct_n_bases'], '31.29', 'pct_n_bases is 31.29')
        self.assertEqual(qc_summary['pct_covered_bases'], '68.01', 'pct_covered_bases is 68.01')
        self.assertEqual(qc_summary['total_variants'], 10, 'total_variants is correct')
        self.assertEqual(qc_summary['total_snv'], 9, 'total_snv is correct')
        self.assertEqual(qc_summary['total_indel'], 1, 'total_indel is correct')
        self.assertEqual(qc_summary['total_n'], 2, 'total_n is correct')
        self.assertEqual(qc_summary['total_iupac'], 4, 'total_iupac is correct')
        self.assertEqual(qc_summary['mean_depth'], 679.4, 'mean_depth is correct')
        self.assertEqual(qc_summary['median_depth'], 682, 'median_depth is correct')
        self.assertEqual(qc_summary['ct'], '17.4', 'ct value is correct')
        self.assertEqual(qc_summary['qc_pass'], 'FALSE', 'qc_pass is correct')

    def test_create_qc_summary_line_no_meta(self):
        qcfile = 'data/sampleA.qc.csv'
        metafile = 'data/metadata.tsv'
        varfile = 'data/sampleA.variants.tsv'
        covfile = 'data/sampleA.per_base_coverage.bed'
        qc_summary = create_qc_summary_line(var_file=varfile,
                                            qc_file=qcfile,
                                            cov_file=covfile,
                                            indel=True)
        self.assertEqual(qc_summary['sample_name'], 'sampleA', 'Sample name correctly identified as: sampleA')
        self.assertEqual(qc_summary['pct_n_bases'], '31.29', 'pct_n_bases is 31.29')
        self.assertEqual(qc_summary['pct_covered_bases'], '68.01', 'pct_covered_bases is 68.01')
        self.assertEqual(qc_summary['total_variants'], 10, 'total_variants is correct')
        self.assertEqual(qc_summary['total_snv'], 9, 'total_snv is correct')
        self.assertEqual(qc_summary['total_indel'], 1, 'total_indel is correct')
        self.assertEqual(qc_summary['total_n'], 2, 'total_n is correct')
        self.assertEqual(qc_summary['total_iupac'], 4, 'total_iupac is correct')
        self.assertEqual(qc_summary['mean_depth'], 679.4, 'mean_depth is correct')
        self.assertEqual(qc_summary['median_depth'], 682, 'median_depth is correct')
        self.assertEqual(qc_summary['ct'], 'NA', 'ct value is correct')
        self.assertEqual(qc_summary['qc_pass'], 'FALSE', 'qc_pass is correct')

    def test_count_iupac_in_fasta(self):
        iupac_file = 'data/tester.fa'
        iupac_count = count_iupac_in_fasta(fasta=iupac_file)
        self.assertEqual(iupac_count, 9, 'IUPAC count is correct')


if __name__ == '__main__':
    unittest.main()
