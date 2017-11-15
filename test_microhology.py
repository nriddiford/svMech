import unittest

from microhom import (
    microhomology,
    longestMatch,
    run_script,
    reversed_seq
)

class MicroHomTestCase(unittest.TestCase):
    """Tests for `microhom.py`."""

    # microhomology()
    def test_is_microhomology(self):
        """Is microhomology successfully detected?"""
        (position, longest_hom, mhseq) =  microhomology('GGGGGTTTTTTTTTTCCCAA', 'CCCAATTTTTTTTTTGGGGG')
        self.assertTrue(longest_hom == 5)
        self.assertTrue(mhseq == 'CCCAA')
        (position, longest_hom, mhseq) =  microhomology('GGGGGTTTTTTTTTTACCCC', 'CCCCATTTTTTTTTTGGGGG')
        self.assertTrue(longest_hom == 4)


    def test_is_no_microhomology(self):
        """Is lack of microhomology successfully reported?"""
        (position, longest_hom, mhseq) =  microhomology('GGGGGTTTTTTTTTTCCCGG', 'CCCAATTTTTTTTTTGGGGG')
        self.assertTrue(longest_hom == 0)


    # reversed_seq()
    def test_reversed_seq(self):
        """Does reversed_seq() successfully reverse these sequences?"""
        self.assertEqual(reversed_seq('ATGC'), 'CGTA')
        self.assertEqual(reversed_seq('AAAA'), 'AAAA')


    # Tests for longestMatch
    def test_homology(self):
        """Does longestMatch() successfully report matches in these sequences?"""
        (downstream_start, downstream_end, upstream_start, upstream_end, homseq) = longestMatch('AAACCCCCCCAAA','GGGCCCCCCCGGG')
        self.assertEqual(len(homseq), 7)

    def test_no_homology(self):
        """Does longestMatch() successfully report matches in these sequences?"""
        (downstream_start, downstream_end, upstream_start, upstream_end, homseq) = longestMatch('AAAAAAAAAA','GGGGGGGGGG')
        self.assertEqual(len(homseq), 0)


    def test_no_hom_no_del_no_ins(self):
        (position, longest_hom, mhseq) =  microhomology('TTCCCC', 'ATGGC')
        self.assertEqual(len(mhseq), 0)

        """
        upstream_seq   = 'GGGAATTTTTTTTTTCCCAT'
        downstream_seq = 'CCCAAGGGTTTTTTTTTTGG'
        """

        (longest_hom, mhseq, homseq, delsize, delbases) = run_script('', '', 'TTTTTCCCATCCCAAGGG', '', 1,0)
        self.assertTrue(delsize == 0)


    def test_microhomology_no_deletion(self):
        """Does script() successfully report matches in these sequences?

        upstream_seq   = 'GGGAATTTTTTTTTTCCCAA'
        downstream_seq = 'CCCAAGGGTTTTTTTTTTGG'

        """
        (longest_hom, mhseq, homseq, delsize, delbases) = run_script('', '', 'TTTTTCCCAAGGGTTTTT', '', 1,1)
        self.assertEqual(longest_hom, 5)
        self.assertTrue(mhseq == 'CCCAA')
        self.assertFalse(mhseq == 'CACAG')
        self.assertTrue(homseq == 'TTTTTTTTTT')
        self.assertEqual(len(homseq), 10)

    def test_microhomology_deletion(self):
        """Does script() successfully report matches in these sequences?

        upstream_seq   = 'GGGAATTTTTTTTTTCCCAA'
        downstream_seq = 'CCCAAGGGTTTTTTTTTTGG'

        """
        (longest_hom, mhseq, homseq, delsize, delbases) = run_script('', '', 'TTTTTCCCAAGGGTTTTT', '', 1,1)
        self.assertTrue(delsize == 0)

        (longest_hom, mhseq, homseq, delsize, delbases) = run_script('', '', 'TTTTTCCCGGGTTTTT', '', 1,1)
        self.assertTrue(delsize == 2)
        self.assertTrue(delbases == 'AA')






if __name__ == '__main__':
    unittest.main()
