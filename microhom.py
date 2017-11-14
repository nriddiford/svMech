#!/usr/bin/python

import pysam
import sys
import re
from difflib import SequenceMatcher

from Bio.Seq import Seq

import argparse

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-p', action="store", help="Breakpoint position (chr:bp1-bp2) [Required]", dest='location', required=True)
    parser.add_argument('-s', action="store", help="Split read sequence for detection of insertions", dest="split")
    parser.add_argument('-n', action="store", help="Number of bases to look for extended homology. [Default: 200 if microhomology found, 10 if not]", dest='homspace',type=int)
    parser.add_argument('-t', action="store_true", help="Run in test mode [Overides -p & -s values]", dest='test')

    args = parser.parse_args()
    pos = args.location
    # homspace = args.homspace
    split_read = args.split
    n=args.homspace
    test=args.test
    genome = pysam.Fastafile("/Users/Nick_curie/Documents/Curie/Data/Genomes/Dmel_v6.12/Dmel_6.12.fasta")
    run_script(pos, n, split_read, genome,test)

def reversed_seq(x):
    """
    Reverses a sequence `x`

    >>> reversed_seq('ATGC') == 'CGTA'
    True
    """
    return x[::-1]

def get_parts(lookup):
    (chrom1, bp1, bp2) = re.split(':|-',lookup)
    return(chrom1, int(bp1), int(bp2))

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
            return(position , longest_hom, mhseq)
        position += 1
        i -= 1
    if longest_hom == 0:
        return False


def getMechanism(homlen, inslen):
    if inslen >= 10:
        mechanism="FoSTeS"
    elif (homlen <= 2) or (inslen >=1 and inslen <= 10):
        mechanism="NHEJ"
    elif homlen >= 3 and homlen <= 100:
        mechanism="Alt-EJ"
    elif homlen >= 100:
        mechanism="NAHR"
    else:
        mechanism="NA"

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

    print("Seq1:%s (%s bases upstream from breakpoint)") % (seq1, len(seq1))
    print("Seq2:%s (%s bases downstream from breakpoint)") % (seq2, len(seq2))

    return(seq1_start, seq1_end, seq2_start, seq2_end, seq)

def run_script(pos, n, split_read, genome,test):
    if not test:
        (chrom1, bp1, bp2) = get_parts(pos)

        # bp1 -= 1
        bp2 -= 1
        upstream = bp1 - 100
        downstream = bp2 + 100

        upstream_seq = genome.fetch(chrom1, upstream, bp1)
        downstream_seq = genome.fetch(chrom1, bp2, downstream)
        # upstream_seq = reversed_seq(upstream_n)
    else:
        # with mh
        upstream_seq = 'GGGGGTTTTTTTTTTCCCAA'
        downstream_seq = 'CCCAATTTTTTTTTTGGGGG'

        # without mh
        # upstream_seq = 'GGGGGTTTTTTTTTTCCCAA'
        # downstream_seq = 'GGGAATTTTTTTTTTGGGGG'

        # split_read = 'TTTTTCCCAAGGGAATTTTT'

        # split_read = reversed_seq(split_read)

    ##################
    ## Microomology ##
    ##################

    print("Upstream:   %s") % (upstream_seq)
    print("Downstream: %s") % (downstream_seq)

    if microhomology(upstream_seq, downstream_seq):
        (position, longest_hom, mhseq) =  microhomology(upstream_seq, downstream_seq)
    else:
        (position, longest_hom, mhseq) = (0,0,'')

    if(longest_hom>0):
        print("")
        print("* Microhomology at breakpoint: %s (%s bp)") % (mhseq, longest_hom)
        downstream_spacer = (len(upstream_seq) - len(mhseq))
        dmarker = " "*(downstream_spacer) + "^"*len(mhseq)
        downstream_buff = " "*downstream_spacer

        # upstream_base_end = bp1 - len(mhseq)
        # downstream_base_end = bp2 + len(mhseq)
        marker = "^"*len(mhseq)
        print(" Upstream:      %s") % (upstream_seq)
        print(" Downstream:    %s%s") % (downstream_buff, downstream_seq)
        print(" Microhomology: %s\n") % (dmarker)

    else:
        print("\n* No microhomology found\n")

    if n is None:
        if longest_hom<=3:
            n=10
        else:
            n=200

    ##############
    ## Homology ##
    ##############

    if n >= len(upstream_seq):
         n = len(upstream_seq)

    # upstream_seq[-n:] takes n characters from the end of the string
    (downstream_start, downstream_end, upstream_start, upstream_end, homseq) = longestMatch(downstream_seq[0:n],upstream_seq[-n:])
    if len(homseq)>0:

        downstream_spacer = len(upstream_seq)  - n  +  upstream_start - downstream_start
        print("* Longest homologous sequence +/- %s bp from breakpoint: %s (%s bp)\n") % (n, homseq, len(homseq))

        seq_diff = len(upstream_seq) - len(upstream_seq[-n:])

        umarker = " "*(upstream_start+seq_diff) + "^"*len(homseq)
        dmarker = " "*(downstream_start+downstream_spacer) + "^"*len(homseq)

        padded_downstream = " "*downstream_spacer + downstream_seq

        print(" Upstream:      %s") % (upstream_seq)
        print(" Homology:      %s") % (umarker)
        print(" Downstream:    %s") % (padded_downstream)
        print(" Homology:      %s\n") % (dmarker)

    else:
        print("\n* No homology found +/- %s bp from breakpoint\n") % (n)


    if split_read is not None:
        ###############
        ## Inserions ##
        ###############

        print("* Split read aligned to upstream sequence:")
        # Align split read to upstream seq (normal orientation)
        (upstream_start, upstream_end, split_start_up, split_end_up, upseq) = longestMatch(upstream_seq, split_read)
        print(upstream_start, upstream_end, split_start_up, split_end_up)
        aligned_up = len(upseq)
        deletion_size= len(upstream_seq[upstream_end:])

        # microhomolgy?
        if longest_hom > 0:

            if deletion_size > 0:
                deleted_bases = "."*(deletion_size)
                print("deletion of %s bases(%s)" % (deletion_size, deleted_bases))

                print(" Upstream:      %s") % (upstream_seq[upstream_start:])
                print(" Split read:    %s%s--/--%s\n") % (split_read[split_start_up:split_end_up], deleted_bases, split_read[split_end_up:len(split_read)])

                (downstream_start, downstream_end, split_start, split_end, downseq) = longestMatch(downstream_seq, split_read)
                print(downstream_start, downstream_end, split_start, split_end, downseq)
                # if longest_hom > 0 and split_read[split_start_up:split_end_up] == mhseq[:-deleted_bases]:

                # difference = downstream_start - split_start
                # print(difference)

                print(mhseq[:-deletion_size])
                add2split = mhseq[:-deletion_size]
                # print(split_read[split_end-add2split:split_end])
                print(split_end_up, longest_hom, deletion_size)
                difference = (split_end_up - longest_hom) + deletion_size
                print(difference)
                print(split_read[difference:])

                seqbuffer = " "*(5+len(split_read[0:split_end_up-len(add2split)]))
                print(" Split read:    %s--/--%s%s%s") % (split_read[0:split_end_up-len(add2split)],deleted_bases, add2split, split_read[split_end_up:])
                print(" Downstream:    %s%s\n") % (seqbuffer, downstream_seq[0:downstream_end])


            else:
                print("\n* No deletion at breakpoint")

                print(" Upstream:      %s") % (upstream_seq[upstream_start:])
                print(" Split read:    %s--/--%s\n") % (split_read[split_start_up:split_end_up], split_read[split_end_up:len(split_read)])


        else:
            if deletion_size > 0:
                deleted_bases = "."*(deletion_size)
                print("deletion of %s bases(%s)" % (deletion_size, deleted_bases))
                print(" Upstream:      %s") % (upstream_seq[upstream_start:])
                print(" Split read:    %s--/--%s\n") % (split_read[split_start_up:split_end_up], split_read[split_end_up:len(split_read)])

            else:
                print("\n* No microhomolgy found, and no deletion or insertion at breakpoint")

                print(" Upstream:      %s") % (upstream_seq[upstream_start:])
                print(" Split read:    %s--/--%s\n") % (split_read[split_start_up:split_end_up], split_read[split_end_up:len(split_read)])

                (downstream_start, downstream_end, split_start, split_end, downseq) = longestMatch(downstream_seq, split_read)
                difference = (downstream_start-split_start)

                seqbuffer = " "*(5+len(split_read[0:split_end])+difference)
                nonaligned = "-"*difference
                print(" Split read:    %s--/--%s%s") % (split_read[0:split_end_up],nonaligned, split_read[split_end_up:])
                print(" Downstream:    %s%s\n") % (seqbuffer, downstream_seq[0:downstream_end])


        bp_start = split_end_up # for extracting the inserted seq


        # Align split read to downstream seq (normal orientation)
        # seq = Seq(split_read)
        ###
        #split_read = seq.reverse_complement()
        insertion_size = 0
        deletion_size = 0

        (downstream_start, downstream_end, split_start, split_end, downseq) = longestMatch(downstream_seq, split_read)
        # if downstream_start == split_start:
        #     print("perfect match")
        #     seqbuffer = " "*(5+len(split_read[0:bp_start-longest_hom]))
        #
        #     print(" Split read:    %s--/--%s%s") % (split_read[0:bp_start-longest_hom],mhseq, split_read[bp_start:])
        #     print(" Downstream:    %s%s\n") % (seqbuffer, downstream_seq[0:downstream_end])
        #
        # else:
        #     print("Not a perfect match")
        #     (downstream_start, downstream_end, split_start, split_end, downseq) = longestMatch(downstream_seq, split_read[bp_start:])
        #     print(downstream_start, downstream_end, split_start, split_end, downseq)
        #
        #     if downstream_start > split_start:
        #         difference = downstream_start - split_start
        #         print(difference)
        #         print(split_read[bp_start-difference:bp_start])
        #
        #         if longest_hom > 0 and mhseq == split_read[bp_start-difference:bp_start]:
        #             print(mhseq)
        #             print("Perfect match with microhomology")
        #             seqbuffer = " "*(5+len(split_read[0:bp_start-longest_hom]))
        #
        #             print(" Split read:    %s--/--%s%s") % (split_read[0:bp_start-longest_hom],mhseq, split_read[bp_start:])
        #             print(" Downstream:    %s%s\n") % (seqbuffer, downstream_seq[0:downstream_end])
        #
        #         else:
        #             print("Deletion scenario")
        #             difference = downstream_start - split_start
        #             print(difference)
        #             deleted_bases = "."*(downstream_start)
        #             deletion_size = len(deleted_bases)
        #             seqbuffer = " "*(5+len(split_read[0:bp_start]))
        #
        #             print(" Split read:    %s--/--%s%s") % (split_read[0:bp_start],deleted_bases, split_read[bp_start:])
        #             print(" Downstream:    %s%s\n") % (seqbuffer, downstream_seq[0:downstream_end])
        #
        #             if deletion_size:
        #                 deleted_seq = downstream_seq[0:deletion_size]
        #                 print("* %s bp deletion '%s' at breakpoint\n") % (deletion_size, deleted_seq)
        #     else:
        #         print("Insertion scenario")

    #     # Calculate length of aligned sequences
    #     split_len = len(split_read)
    #     aligned_up = len(upseq)
    #     aligned_down = len(downseq)
    #     print(aligned_up, aligned_down, split_len)
    #     print(downstream_start, downstream_end, split_start, split_end)
    #     # Split read length - aligned portion = insetion size
    #     insertion_size = int(split_len - aligned_up - aligned_down)
    #
    #     print(insertion_size)
    #     inserted_seq = split_read[bp_start:bp_start+insertion_size]
    #
    #     print("* Split read aligned to downstream sequence:")
    #     if insertion_size <= 0:
    #         deleted_bases = "."*(downstream_start)
    #         deletion_size = len(deleted_bases)
    #         seqbuffer = " "*(5+len(split_read[0:bp_start]))
    #
    #         print(" Split read:    %s--/--%s%s") % (split_read[0:bp_start],deleted_bases, split_read[bp_start:])
    #         print(" Downstream:    %s%s\n") % (seqbuffer, downstream_seq[0:downstream_end])
    #
    #         if deletion_size:
    #             deleted_seq = downstream_seq[0:deletion_size]
    #             print("* %s bp deletion '%s' at breakpoint\n") % (deletion_size, deleted_seq)
    #
    #     else:
    #         seqbuffer = " "*(downstream_start+5+insertion_size)
    #
    #         print(" Split read:    %s--/--%s") % (split_read[0:split_end-1], split_read[split_end-1:])
    #         print(" Downstream:    %s     %s\n") % (seqbuffer, downstream_seq[downstream_start:])
    #
    #         print("* %s bp insertion '%s' at breakpoint\n") % (insertion_size, inserted_seq)
    #
    #     if insertion_size >= 3:
    #         (inserted_start, inserted_end, upstream_start, upstream_end, aligned) = longestMatch(inserted_seq,upstream_seq)
    #         if len(aligned) > 3:
    #             splitbuffer = " "*upstream_start
    #             print(" Upstream:      %s--/--") % (upstream_seq)
    #             print(" Insertion:     %s%s\n") % (splitbuffer, inserted_seq[0:len(aligned)])
    #             insertion_pos = (len(upstream_seq) - len(aligned))
    #             print("* %s bp of inserted sequence -%s bps from breakpoint on upstream sequence\n") % (len(aligned),insertion_pos)
    #
    #         (downstream_start, downstream_end, inserted_start, inserted_end, aligned) = longestMatch(downstream_seq, inserted_seq)
    #         if len(aligned) > 3:
    #             print(downstream_start)
    #             splitbuffer = " "*(downstream_start+5)
    #             print(" Downstream:     --/--%s") % (downstream_seq)
    #             print(" Insertion:      %s%s\n") % (splitbuffer, inserted_seq)
    #             print("* %s bp of inserted sequence +%s bps from breakpoint on downstream sequence\n") % (len(aligned),downstream_start)
    #
    # else:
    #     print("No split read sequence provided. Unable to find insetion at bp")
    #     insertion_size = 0
    #
    # ###############
    # ## Mechanism ##
    # ###############
    #
    # # Calculate mechanism
    # mechanism=getMechanism(longest_hom, insertion_size)
    # print("* Mechanism: %s") % (mechanism)
    # if deletion_size >= 0:
    #     print("  * %s bp deletion") % (deletion_size)
    # else:
    #     print("  * %s bp insertion") % (insertion_size)
    # print("  * %s bp homology at breakpoints") % (longest_hom)
    #
    #
    #
    # #
    return(longest_hom, mhseq, homseq)


if __name__ == '__main__':
    main()
