import unittest

from microhom import microhomology
from microhom import reversed_seq

class MicroHomTestCase(unittest.TestCase):
    """Tests for `microhom.py`."""

    # microhomology()
    def test_is_microhomology(self):
        """Is microhomology successfully detected?"""
        self.assertTrue(microhomology('GGGGGTTTTTTTTTTCCCAA', 'CCCAATTTTTTTTTTGGGGG'))
        self.assertFalse(microhomology('GGGGGTTTTTTTTTTCCCGG', 'CCCAATTTTTTTTTTGGGGG'))

    def test_is_no_microhomology(self):
        """Is lack of microhomology successfully reported?"""
        self.assertFalse(microhomology('GGGGGTTTTTTTTTTCCCGG', 'CCCAATTTTTTTTTTGGGGG'))
        # (position, longest_hom, mhseq) =  microhomology('GGGGGTTTTTTTTTTCCCGG', 'CCCAATTTTTTTTTTGGGGG')
        # self.assertEqual(longest_hom, False)

    def test_microhomology_length(self):
        """Is length of microhomology successfully reported?"""
        (position, longest_hom, mhseq) =  microhomology('GGGGGTTTTTTTTTTCCCAA', 'CCCAATTTTTTTTTTGGGGG')
        self.assertEqual(longest_hom, 5)
        (position, longest_hom, mhseq) =  microhomology('GGGGGTTTTTTTTTT', 'TTTTTTTTGGGGG')
        self.assertEqual(longest_hom, 8)

    # reversed_seq()
    def test_reversed_seq(self):
        """Does reversed_seq() successfully reverse these sequences?"""
        self.assertEqual(reversed_seq('ATGC'), 'CGTA')
        self.assertEqual(reversed_seq('AAAA'), 'AAAA')

if __name__ == '__main__':
    unittest.main()
