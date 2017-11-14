import unittest

from microhom import microhomology, longestMatch, run_script
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


    # Tests for longestMatch
    def test_match(self):
        """Does longestMatch() successfully report matches in these sequences?"""
        (downstream_start, downstream_end, upstream_start, upstream_end, seq) = longestMatch('AAACCCCCCCAAA','GGGCCCCCCCGGG')
        self.assertEqual(len(seq), 7)

    def test_no_match(self):
        """Does longestMatch() successfully report matches in these sequences?"""
        (downstream_start, downstream_end, upstream_start, upstream_end, seq) = longestMatch('AAAAAAAAAA','GGGGGGGGGG')
        self.assertEqual(len(seq), 0)



    def test_script(self):
        """Does longestMatch() successfully report matches in these sequences?"""
        (longest_hom, mhseq, homseq) = run_script('X:1-2', 20, 'TTTTTCCCCCTTTTT', genome='',test=1)
        self.assertEqual(longest_hom, 5)
        self.assertTrue(mhseq == 'CCCAA')
        self.assertTrue(homseq == 'TTTTTTTTTT')
        self.assertEqual(len(homseq), 10)




        # self.assertEqual(len(seq), 0)


if __name__ == '__main__':
    unittest.main()
