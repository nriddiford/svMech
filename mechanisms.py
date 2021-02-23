#!/usr/bin/python

import pysam
import sys
import re
from optparse import OptionParser
from difflib import SequenceMatcher
from string import maketrans

from Bio.Seq import Seq


def get_args():
    parser = OptionParser()
    parser.add_option('-g', action="store", help="Genome fasta file", dest='genomeFile')
    parser.add_option('-p', action="store", help="Breakpoint position (chr:bp1-bp2) [Required]", dest='location')
    parser.add_option('-s', action="store", help="Split read sequence for detection of insertions", dest="split")
    parser.add_option('-n', action="store", help="Number of bases to look for extended homology. [Default: 200 if microhomology found, 10 if not]", dest='homspace',type=int)
    parser.add_option('-o', action="store", help="Orientation of split read [Default FF. Options FR RF]", dest="ori")
    parser.add_option('-c', action="store", help="Amount of upstream/downstream seq to consider [Default full length]", dest="cut")

    parser.set_defaults(ori = 'FF',
                        cut = 100,
                        homspace = 200,
                        genomeFile = "/Users/Nick_curie/Documents/Curie/Data/Genomes/Dmel_v6.12/Dmel_6.12.fasta")

    options, args = parser.parse_args()

    if not options.location:
        parser.print_help()
        print
    return (options, args)


def reverse_complement(seq):
    return seq.translate(maketrans('ACGTacgt', 'TGCAtgca'))


def reversed_seq(x):
    return x[::-1]


def get_parts(lookup):
    (chrom, bp1, bp2) = re.split(':|-',lookup)
    return(chrom, int(bp1), int(bp2))


def microhomology(seq1, seq2):
    mh = 0
    pos_count = 0
    mhseq = ""
    position =""
    longest_hom=0
    i = len(seq1)
    position = 0
    while i > 0:
        # print("Up   %s: %s") % (i, seq1[position:])
        # print("Down %s: %s") % (i, seq2[:i])
        if seq1[position:] == seq2[:i]:
            position = i
            longest_hom = len(seq2[:i])
            mhseq = seq2[:i]
        position += 1
        i -= 1

    if(longest_hom>0):
        return(position, longest_hom, mhseq)
    else:
        return(0, 0, '')


def getMechanism(homlen, inslen, delsize, templen):
    if homlen >= 100:
        mechanism="NAHR"
    elif inslen >= 10 or templen >= 5:
        mechanism="FoSTeS"
    elif homlen >= 3 and homlen <= 100 and inslen == 0:
        mechanism="Alt-EJ"
    else:
        mechanism="NHEJ"

    return(mechanism)


def longestMatch(seq1, seq2):
    s = SequenceMatcher(None, seq1, seq2)
    match = s.find_longest_match(0, len(seq1), 0, len(seq2))

    block = s.get_matching_blocks()
    seq1_start = match[0]
    seq1_end = match[0]+match[2]
    seq = seq1[match[0]:(match[0]+match[2])]
    seq2_start = match[1]
    seq2_end = match[1]+match[2]
    # print("Seq1:%s") % (seq1)
    # print("Seq2:%s") % (seq2)
    # print(seq1_start, seq1_end, seq2_start, seq2_end, seq)
    return(seq1_start, seq1_end, seq2_start, seq2_end, seq)

def getSeqs(chrom, bp1, bp2, homologySearch, cut, genome):
    bp2 -= 1
    upstram_start = bp1 - cut
    upstream_of_bp1 = bp1 - 50
    downstream_of_bp1 = bp1 + 50
    upstream_seq = genome.fetch(chrom, upstram_start, bp1)
    bp1_surrounding = genome.fetch(chrom, upstream_of_bp1, downstream_of_bp1)

    print("Upstream:   %s") % (upstream_seq)

    downstream_stop = bp2 + cut
    upstream_of_bp2 = bp2 - 50
    downstream_of_bp2 = bp2 + 50
    downstream_seq = genome.fetch(chrom, bp2, downstream_stop)
    bp2_surrounding = genome.fetch(chrom, upstream_of_bp2, downstream_of_bp2)

    print("Downstream: %s") % (downstream_seq)

    if (bp2-bp1<=homologySearch):
        homologySearch = bp2-bp1

    inside_window = bp2 - homologySearch
    inside_seq = genome.fetch(chrom, inside_window, bp2)
    print("Inside: %s") % (inside_seq)

    return inside_seq, upstream_seq, downstream_seq




def main():
    options, args = get_args()
    genome = pysam.Fastafile(options.genomeFile)

    pos = options.location
    split_read = options.split
    homologySearch = options.homspace
    ori = options.ori
    cut = options.cut

    (chrom, bp1, bp2) = get_parts(pos)

    iSeq, uSeq, dSeq = getSeqs(chrom, bp1, bp2, homologySearch, cut, genome)

    position, longest_hom, mhseq = microhomology(uSeq, iSeq)

    if longest_hom > 0:
        print("\n* Microhomology at breakpoint: %s (%s bp)") % (mhseq, longest_hom)
    #     downstream_spacer = (len(upstream_seq[-50:]) - len(mhseq))
    #     dmarker = " "*(downstream_spacer) + "^"*len(mhseq)
    #     downstream_buff = " "*downstream_spacer
    #
    #     # upstream_base_end = bp1 - len(mhseq)
    #     # downstream_base_end = bp2 + len(mhseq)
    #     marker = "^"*len(mhseq)
    #     print(" Upstream:      %s") % (upstream_seq[-50:])
    #     print(" Downstream:    %s%s") % (downstream_buff, downstream_seq[:50])
    #     print(" Microhomology: %s\n") % (dmarker)
    #




if __name__ == "__main__":
    sys.exit(main())
