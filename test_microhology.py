import unittest

from microhom import (
    microhomology,
    longestMatch,
    run_script,
    reverse_complement
)

class MicroHomTestCase(unittest.TestCase):

    def test_is_microhomology(self):
        """ Is microhomology successfully detected? """
        (position, longest_hom, mhseq) =  microhomology('GGGGGTTTTTTTTTTCCCAA', 'CCCAATTTTTTTTTTGGGGG')
        self.assertTrue(longest_hom == 5)
        self.assertTrue(mhseq == 'CCCAA')
        (position, longest_hom, mhseq) =  microhomology('GGGGGTTTTTTTTTTACCCC', 'CCCCATTTTTTTTTTGGGGG')
        self.assertTrue(longest_hom == 4)

    def test_microhomology_no_deletion(self):
        """Is the microhomology correctly reported in the sequences? """

        upstream_seq   = 'GGGAATTTTTTTTTTCCCAA'
        downstream_seq = 'CCCAAGGGTTTTTTTTTTGG'
        split_read     = 'TTTTTCCCAAGGGTTTTT'

        (longest_hom, mhseq, homseq, delsize, delbases, insize, inseq, tempup, tempdown, templen, mech, seq1, seq2) = run_script('', '', split_read, '', upstream_seq, downstream_seq, 'FF', 200)
        self.assertEqual(longest_hom, 5)
        self.assertTrue(mhseq == 'CCCAA')
        self.assertFalse(mhseq == 'CACAG')
        self.assertTrue(homseq == 'TTTTTTTTTT')
        self.assertTrue(len(homseq) == 10)


class HomTestCase(unittest.TestCase):

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

class NoMicroHomTestCase(unittest.TestCase):

    def test_is_no_microhomology(self):
        """ Is lack of microhomology successfully reported? """
        (position, longest_hom, mhseq) =  microhomology('GGGGGTTTTTTTTTTCCCGG', 'CCCAATTTTTTTTTTGGGGG')
        self.assertTrue(longest_hom == 0)

    def test_no_hom_no_del_no_ins(self):
        """ Are these perfect matches with no mh reported correctly? """

        (position, longest_hom, mhseq) =  microhomology('TTCCCC', 'ATGGCT')
        self.assertEqual(len(mhseq), 0)

        upstream_seq   = 'TTTTTTT'
        downstream_seq = 'CCCCCCC'
        split_read     = 'TTTTTTTCCCCCCC'

        (longest_hom, mhseq, homseq, delsize, delbases, insize, inseq, tempup, tempdown, templen, mech, seq1, seq2) = run_script('', '', split_read, '', upstream_seq, downstream_seq, 'FF', 200)
        self.assertTrue(delsize == 0)
        self.assertTrue(insize == 0)

        upstream_seq   = 'GGGAATTTTTTTTTTCCCAT'
        downstream_seq = 'CCCAAGGGTTTTTTTTTTGG'
        split_read     = 'TTTTTCCCATCCCAAGGG'

        (longest_hom, mhseq, homseq, delsize, delbases, insize, inseq, tempup, tempdown, templen, mech, seq1, seq2) = run_script('', '', split_read, '', upstream_seq, downstream_seq, 'FF', 200)
        self.assertTrue(delsize == 0)
        self.assertTrue(insize == 0)


class DeletionTestCase(unittest.TestCase):

    def test_microhomology_deletion(self):
        """Is the deletion with microhomology correctly reported in these sequences? """

        upstream_seq   = 'GGGAATTTTTTTTTTCCCAA'
        downstream_seq = 'CCCAAGGGTTTTTTTTTTGG'
        split_read     = 'TTTTTCCCAAGGGTTTTT'

        (longest_hom, mhseq, homseq, delsize, delbases, insize, inseq, tempup, tempdown, templen, mech, seq1, seq2) = run_script('', '', split_read, '', upstream_seq, downstream_seq, 'FF', 200)
        self.assertTrue(delsize == 0)

        split_read     = 'TTTTTCCCGGGTTTTT'

        (longest_hom, mhseq, homseq, delsize, delbases, insize, inseq, tempup, tempdown, templen, mech, seq1, seq2) = run_script('', '', split_read, '', upstream_seq, downstream_seq, 'FF', 200)
        self.assertTrue(delsize == 2)
        self.assertTrue(delbases == 'AA')

    def test_no_microhomology_deletion(self):
        """Is the deletion with no microhomology correctly reported in these sequences? """

        upstream_seq   = 'AAAATTTT'
        downstream_seq = 'GGGGCCCC'
        split_read     = 'AAAATTTTGGCCCC'

        (longest_hom, mhseq, homseq, delsize, delbases, insize, inseq, tempup, tempdown, templen, mech, seq1, seq2) = run_script('', '', split_read, '', upstream_seq, downstream_seq, 'FF', 200)
        self.assertTrue(delsize == 2)
        self.assertTrue(delbases == 'GG')

        upstream_seq   = 'ATATGGCC'
        downstream_seq = 'GCATCCCTTT'
        split_read     = 'ATATGGGCATCCCTTT'

        (longest_hom, mhseq, homseq, delsize, delbases, insize, inseq, tempup, tempdown, templen, mech, seq1, seq2) = run_script('', '', split_read, '', upstream_seq, downstream_seq, 'FF', 200)
        self.assertTrue(delsize == 2)
        self.assertTrue(delbases == 'CC')


class InsertionTestCase(unittest.TestCase):

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

        (longest_hom, mhseq, homseq, delsize, delbases, insize, inseq, tempup, tempdown, templen, mech, seq1, seq2) = run_script('', '', split_read, '', upstream_seq, downstream_seq, 'FF', 200)
        self.assertTrue(insize == 6)
        self.assertTrue(inseq == 'TAGCAT')

        split_read = 'TTTTTCCCAATACCCAAGGGTTTTT'
        (longest_hom, mhseq, homseq, delsize, delbases, insize, inseq, tempup, tempdown, templen, mech, seq1, seq2) = run_script('', '', split_read, '', upstream_seq, downstream_seq, 'FF', 200)
        self.assertTrue(insize == 2)
        self.assertTrue(inseq == 'TA')

        upstream_seq   = 'GGGAATTTTTTAATG'
        downstream_seq = 'ATGGGGTTTTTTGG'
        split_read = 'TTTTAATGCCCATGGGGTT'

        (longest_hom, mhseq, homseq, delsize, delbases, insize, inseq, tempup, tempdown, templen, mech, seq1, seq2) = run_script('', '', split_read, '', upstream_seq, downstream_seq, 'FF', 200)
        self.assertTrue(insize == 3)
        self.assertTrue(inseq == 'CCC')

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

        (longest_hom, mhseq, homseq, delsize, delbases, insize, inseq, tempup, tempdown, templen, mech, seq1, seq2) = run_script('', '', split_read, '', upstream_seq, downstream_seq, 'FF', 200)
        self.assertTrue(tempup == 'TAG')
        self.assertTrue(tempdown == 'GCAT')

    def test_microhomology_templated_insertion2(self):
        """Are the templated insertions in these reads correctly reported?"""
        upstream_seq   = 'GATCCCCCG'
        downstream_seq = 'AAAAAGTCCC'
        split_read     = 'GATCCCCCGGATCAAAAAGTCCC'

        (longest_hom, mhseq, homseq, delsize, delbases, insize, inseq, tempup, tempdown, templen, mech, seq1, seq2) = run_script('', '', split_read, '', upstream_seq, downstream_seq, 'FF', 200)
        self.assertTrue(tempup == 'GATC')
        self.assertTrue(insize == 4)

    def test_microhomology_templated_insertion3(self):
        """Are the templated insertions in these reads correctly reported?"""
        upstream_seq   = 'TGTAATACCATTTACCCTATTGATTTGCTGTTCTTATTGTTATATCATTT'
        downstream_seq = 'ATTACATGAAGTGCAGACTTCTCTTAAAAAAAAAAATCCACTTTAATTTT'
        split_read     = 'ACCATTTACCCTATTGATTTGCTGTTCTTATTGTTATATCATTTTATATCATTGTTTCATGATGGATTACATGAAGTGCAGACTTCTCTTAAAAAAAAAAC'

        """ insertion of TATATCATTGTTTCATGATGG """
        (longest_hom, mhseq, homseq, delsize, delbases, insize, inseq, tempup, tempdown, templen, mech, seq1, seq2) = run_script('', '', split_read, '', upstream_seq, downstream_seq, 'FF', 200)
        self.assertTrue(inseq == 'TATATCATTGTTTCATGATGG')
        self.assertTrue(tempup == 'TATATCATT')
        self.assertTrue(tempdown == 'CATGA')
        self.assertTrue(insize == 21)

    def test_microhomology_templated_insertion4(self):
        upstream_seq   = 'GGGAATAGTTTTTTTTTTATGCC'
        downstream_seq = 'ATGCCGGGTTTTTTTTTTGGCAT'
        split_read = 'TTTTTATGCCTTTGGATGCCGGGT'

        (longest_hom, mhseq, homseq, delsize, delbases, insize, inseq, tempup, tempdown, templen, mech, seq1, seq2) = run_script('', '', split_read, '', upstream_seq, downstream_seq, 'FF', 200)
        self.assertTrue(inseq == 'TTTGG')
        self.assertTrue(tempdown == 'TTTGG')
        self.assertTrue(insize == 5)

# class Inversion_FR(unittest.TestCase):
#     def test_microhomology_inversion(self):
#         """These reads are aligned FR - are they correctly aligned? """
#
#         upstream_seq   = 'TAGCCA'
#         downstream_seq = 'ACCGTCC'
#         split_read     = 'TAGCCAGGACGGT'
#
#         (longest_hom, mhseq, homseq, delsize, delbases, insize, inseq, tempup, tempdown, templen, mech,seq1, seq2) = run_script('', '', split_read, '', upstream_seq, downstream_seq, 'FR', 200)
#         self.assertEqual(longest_hom, 1)
#         self.assertTrue(mhseq == 'A')
#
#     def test_reverse_complement(self):
#         """ Are these sequences correctly reversed? """
#         self.assertEqual(reverse_complement('AGGACGGT'), 'ACCGTCCT')
#
#
#     def test_no_microhomology_inversion(self):
#         """These reads are aligned FR - are they correctly aligned? """
#
#         upstream_seq   = 'TAGCCC'
#         downstream_seq = 'ACCGTCC'
#         split_read     = 'TAGCCCGGACGGT'
#
#         (longest_hom, mhseq, homseq, delsize, delbases, insize, inseq, tempup, tempdown, templen, mech,seq1, seq2) = run_script('', '', split_read, '', upstream_seq, downstream_seq, 'FR', 200)
#         self.assertEqual(longest_hom, 0)
#         self.assertFalse(mhseq == 'A')
#
#
# class Inversion_FR_Del(unittest.TestCase):
#     def test_microhomology_inversion_del(self):
#         """These reads are aligned FR, and have a deletion at the breakpoint. Are they correctly aligned? """
#
#         upstream_seq   = 'TAGCCA'
#         downstream_seq = 'ACCGTCC'
#         split_read     = 'TAGCCAGGACG'
#
#         (longest_hom, mhseq, homseq, delsize, delbases, insize, inseq, tempup, tempdown, templen, mech,seq1, seq2) = run_script('', '', split_read, '', upstream_seq, downstream_seq, 'FR', 200)
#         self.assertEqual(delsize, 1)
#         self.assertTrue(mhseq == 'A')
#
#
#         upstream_seq   = 'GACCGCTCTCCATTGGGCAGCGGCGGTTAACAATACCGAAGCGGTGAACA'
#         downstream_seq = 'AAAAAAAAAAAAAGCGGAGAACAAACTCCCGGGCCCTGTGCACTTAATTT'
#         split_read     = 'TTAACAATACCGAAGCGGTGAACATGTTCTCCGCTTTTTTTTTTTTT'
#
#         (longest_hom, mhseq, homseq, delsize, delbases, insize, inseq, tempup, tempdown, templen, mech,seq1, seq2) = run_script('', '', split_read, '', upstream_seq, downstream_seq, 'FR', 200)
#         self.assertTrue(delsize == 0)
#         self.assertTrue(mhseq == 'A')
#
# class Inversion_FR_Ins(unittest.TestCase):
#     def test_microhomology_inversion_del(self):
#         """These reads are aligned FR, and have an insertion at the breakpoint. Are they correctly aligned? """
#
#         upstream_seq   = 'TAGCCA'
#         downstream_seq = 'ACCGTCC'
#         split_read     = 'TAGCCAGGACGGTTCCAC'
#
#         (longest_hom, mhseq, homseq, delsize, delbases, insize, inseq, tempup, tempdown, templen, mech,seq1, seq2) = run_script('', '', split_read, '', upstream_seq, downstream_seq, 'FR', 200)
#         self.assertEqual(insize, 5)
#         self.assertTrue(inseq == 'GTGGA')
#

if __name__ == '__main__':
    unittest.main()
