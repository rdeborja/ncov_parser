import unittest
from ncov.parser.qc import is_indel

class TestIndel(unittest.TestCase):
    def test_is_indel_fail(self):
        self.assertEqual(is_indel("A"), False, "Should be False, A is not an indel")
    def test_is_indel_succeed(self):
        self.assertEqual(is_indel("+A"), True, "Should be True, +A is an indel")

if __name__ == '__main__':
    unittest.main()
