import unittest
from ncov.parser.qc import is_indel, get_total_variants

class TestIndel(unittest.TestCase):
    def test_is_indel_fail(self):
        self.assertEqual(is_indel("A"), False, "Should be False, A is not an indel")
    def test_is_indel_succeed(self):
        self.assertEqual(is_indel("+A"), True, "Should be True, +A is an indel")
    def test_count_variants_no_indel(self):
        self.assertEqual(get_total_variants(file="../../data/sampleA.variants.tsv", indel=False)["total_variants"], 9, "9 single nucleotide variants")
    def test_count_variants_with_indel(self):
        self.assertEqual(get_total_variants(file="../../data/sampleA.variants.tsv", indel=True)["total_variants"], 10, "10 variants")

if __name__ == '__main__':
    unittest.main()
