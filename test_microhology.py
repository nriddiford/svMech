import unittest

from microhom import (
    microhomology,
    longestMatch,
    run_script,
    reversed_seq
)

class MicroHomTestCase(unittest.TestCase):
    # microhomology()
    def test_is_microhomology(self):
        """ Is microhomology successfully detected? """
        (position, longest_hom, mhseq) =  microhomology('GGGGGTTTTTTTTTTCCCAA', 'CCCAATTTTTTTTTTGGGGG')
        self.assertTrue(longest_hom == 5)
        self.assertTrue(mhseq == 'CCCAA')
        (position, longest_hom, mhseq) =  microhomology('GGGGGTTTTTTTTTTACCCC', 'CCCCATTTTTTTTTTGGGGG')
        self.assertTrue(longest_hom == 4)

    def test_is_no_microhomology(self):
        """ Is lack of microhomology successfully reported? """
        (position, longest_hom, mhseq) =  microhomology('GGGGGTTTTTTTTTTCCCGG', 'CCCAATTTTTTTTTTGGGGG')
        self.assertTrue(longest_hom == 0)

    def test_reversed_seq(self):
        """ Are these sequences correctly reversed? """
        self.assertEqual(reversed_seq('ATGC'), 'CGTA')
        self.assertEqual(reversed_seq('AAAA'), 'AAAA')

    # Tests for longestMatch
    def test_homology(self):
        """ Is the extended homology correctly reported in these sequences? """
        seq1 = 'AAACCCCCCCAAA'
        seq2 = 'GGGCCCCCCCGGG'

        (downstream_start, downstream_end, upstream_start, upstream_end, homseq) = longestMatch(seq1,seq2)
        self.assertEqual(len(homseq), 7)

    def test_no_homology(self):
        """ Is the lack of extended homology correctly reported in these sequences? """
        seq1 = 'AAAAAAAAAA'
        seq2 = 'GGGGGGGGGG'

        (downstream_start, downstream_end, upstream_start, upstream_end, homseq) = longestMatch(seq1,seq2)
        self.assertEqual(len(homseq), 0)


    def test_no_hom_no_del_no_ins(self):
        """ Are these perfect matches with no mh reported correctly? """

        (position, longest_hom, mhseq) =  microhomology('TTCCCC', 'ATGGCT')
        self.assertEqual(len(mhseq), 0)

        upstream_seq   = 'TTTTTTT'
        downstream_seq = 'CCCCCCC'
        split_read     = 'TTTTTTTCCCCCCC'

        (longest_hom, mhseq, homseq, delsize, delbases, insize, inseq, tempup, tempdown) = run_script('', '', split_read, '', upstream_seq, downstream_seq)
        self.assertTrue(delsize == 0)
        self.assertTrue(insize == 0)

        upstream_seq   = 'GGGAATTTTTTTTTTCCCAT'
        downstream_seq = 'CCCAAGGGTTTTTTTTTTGG'
        split_read     = 'TTTTTCCCATCCCAAGGG'

        (longest_hom, mhseq, homseq, delsize, delbases, insize, inseq, tempup, tempdown) = run_script('', '', split_read, '', upstream_seq, downstream_seq)
        self.assertTrue(delsize == 0)
        self.assertTrue(insize == 0)


    def test_microhomology_no_deletion(self):
        """Is the microhomology correctly reported in the sequences? """

        upstream_seq   = 'GGGAATTTTTTTTTTCCCAA'
        downstream_seq = 'CCCAAGGGTTTTTTTTTTGG'
        split_read     = 'TTTTTCCCAAGGGTTTTT'

        (longest_hom, mhseq, homseq, delsize, delbases, insize, inseq, tempup, tempdown) = run_script('', '', split_read, '', upstream_seq, downstream_seq)
        self.assertEqual(longest_hom, 5)
        self.assertTrue(mhseq == 'CCCAA')
        self.assertFalse(mhseq == 'CACAG')
        self.assertTrue(homseq == 'TTTTTTTTTT')
        self.assertTrue(len(homseq) == 10)

    def test_microhomology_deletion(self):
        """Is the deletion with microhomology correctly reported in the sequences? """

        upstream_seq   = 'GGGAATTTTTTTTTTCCCAA'
        downstream_seq = 'CCCAAGGGTTTTTTTTTTGG'
        split_read     = 'TTTTTCCCAAGGGTTTTT'

        (longest_hom, mhseq, homseq, delsize, delbases, insize, inseq, tempup, tempdown) = run_script('', '', split_read, '', upstream_seq, downstream_seq)
        self.assertTrue(delsize == 0)

        split_read     = 'TTTTTCCCGGGTTTTT'

        (longest_hom, mhseq, homseq, delsize, delbases, insize, inseq, tempup, tempdown) = run_script('', '', split_read, '', upstream_seq, downstream_seq)
        self.assertTrue(delsize == 2)
        self.assertTrue(delbases == 'AA')

    def test_microhomology_insertion(self):
        """Is the insertion with microhomology correctly reported in the sequences?

        upstream_seq   = 'GGGAATTTTTTTTTTCCCAA'
        downstream_seq = 'CCCAAGGGTTTTTTTTTTGG'
        split_read = 'TTTTTCCCAATAGCATCCCAAGGGTTTTT'

        correct configuration:
        GGGAATTTTTTTTTTCCCAA
                  TTTTTCCCAA---[TAGCAT]---CCCAAGGGTTTTT
                                          CCCAAGGGTTTTTTTTTTGG
                       ^^^^^              ^^^^^
        6 bp insertion (TAGCAT)
        """

        upstream_seq   = 'GGGAATTTTTTTTTTCCCAA'
        downstream_seq = 'CCCAAGGGTTTTTTTTTTGG'
        split_read = 'TTTTTCCCAATAGCATCCCAAGGGTTTTT'

        (longest_hom, mhseq, homseq, delsize, delbases, insize, inseq, tempup, tempdown) = run_script('', '', split_read, '', upstream_seq, downstream_seq)
        self.assertTrue(insize == 6)
        self.assertTrue(inseq == 'TAGCAT')

        split_read = 'TTTTTCCCAATACCCAAGGGTTTTT'
        (longest_hom, mhseq, homseq, delsize, delbases, insize, inseq, tempup, tempdown) = run_script('', '', split_read, '', upstream_seq, downstream_seq)
        self.assertTrue(insize == 2)
        self.assertTrue(inseq == 'TA')


    def test_microhomology_templated_insertion1(self):
        """Are the templated insertions in these split reads correctly reported?

        upstream_seq   = 'GGGAATAGTTTTTTTTTTCCCAA'
        downstream_seq = 'CCCAAGGGTTTTTTTTTTGGCAT'
        split_read = 'TTTTTCCCAATAGCATCCCAAGGGTTTTT'

        correct configuration:
        GGGAATTTTTTTTTTCCCAA
                  TTTTTCCCAA---[TAGCAT]---CCCAAGGGTTTTT
                                          CCCAAGGGTTTTTTTTTTGG
                       ^^^^^              ^^^^^
        6 bp insertion (TAGCAT)
        """

        upstream_seq   = 'GGGAATAGTTTTTTTTTTCCCAA'
        downstream_seq = 'CCCAAGGGTTTTTTTTTTGGCAT'
        split_read = 'TTTTTCCCAATAGCATCCCAAGGGTTTTT'

        (longest_hom, mhseq, homseq, delsize, delbases, insize, inseq, tempup, tempdown) = run_script('', '', split_read, '', upstream_seq, downstream_seq)
        self.assertTrue(tempup == 'TAG')
        self.assertTrue(tempdown == 'GCAT')

    def test_microhomology_templated_insertion2(self):
        """Are the templated insertions in these reads correctly reported?"""
        upstream_seq   = 'GATCCCCCG'
        downstream_seq = 'AAAAAGTCCC'
        split_read     = 'GATCCCCCGGATCAAAAAGTCCC'

        (longest_hom, mhseq, homseq, delsize, delbases, insize, inseq, tempup, tempdown) = run_script('', '', split_read, '', upstream_seq, downstream_seq)
        self.assertTrue(tempup == 'GATC')
        self.assertTrue(insize == 4)

    def test_microhomology_templated_insertion3(self):
        """Are the templated insertions in these reads correctly reported?"""
        upstream_seq   = 'TGTAATACCATTTACCCTATTGATTTGCTGTTCTTATTGTTATATCATTT'
        downstream_seq = 'ATTACATGAAGTGCAGACTTCTCTTAAAAAAAAAAATCCACTTTAATTTT'
        split_read     = 'ACCATTTACCCTATTGATTTGCTGTTCTTATTGTTATATCATTTTATATCATTGTTTCATGATGGATTACATGAAGTGCAGACTTCTCTTAAAAAAAAAAC'

        """ insertion of TATATCATTGTTTCATGATGG """
        (longest_hom, mhseq, homseq, delsize, delbases, insize, inseq, tempup, tempdown) = run_script('', '', split_read, '', upstream_seq, downstream_seq)
        self.assertTrue(inseq == 'TATATCATTGTTTCATGATGG')
        self.assertTrue(tempup == 'TATATCATT')
        self.assertTrue(tempdown == 'CATGA')
        self.assertTrue(insize == 21)



    def test_no_microhomology_deletion(self):
        """Does script() successfully report matches in these sequences?"""

        upstream_seq   = 'GGGAATTTTTTTTTTCCCAA'
        downstream_seq = 'CCCAAGGGTTTTTTTTTTGG'
        split_read     = 'TTTTTCCCAAGGGTTTTT'

        (longest_hom, mhseq, homseq, delsize, delbases, insize, inseq, tempup, tempdown) = run_script('', '', split_read, '', upstream_seq, downstream_seq)
        self.assertTrue(delsize == 0)

        split_read     = 'TTTTTCCCGGGTTTTT'

        (longest_hom, mhseq, homseq, delsize, delbases, insize, inseq, tempup, tempdown) = run_script('', '', split_read, '', upstream_seq, downstream_seq)
        self.assertTrue(delsize == 2)
        self.assertTrue(delbases == 'AA')



if __name__ == '__main__':
    unittest.main()
