#!/usr/bin/python

import pysam
import sys
import re
from difflib import SequenceMatcher
from string import maketrans

from Bio.Seq import Seq

import argparse

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-p', action="store", help="Breakpoint position (chr:bp1-bp2) [Required]", dest='location')
    parser.add_argument('-s', action="store", help="Split read sequence for detection of insertions", dest="split")
    parser.add_argument('-n', action="store", help="Number of bases to look for extended homology. [Default: 200 if microhomology found, 10 if not]", dest='homspace',type=int)
    parser.add_argument('-up', action="store", help="Upstream sequence", dest="upseq")
    parser.add_argument('-down', action="store", help="Upstream sequence", dest="downseq")
    parser.add_argument('-o', action="store", help="Orientation of split read [Default FF. Options FR RF]", dest="ori")

    args = parser.parse_args()
    pos = args.location
    split_read = args.split
    n=args.homspace
    upstream_seq = args.upseq
    downstream_seq = args.downseq
    ori = args.ori
    if ori is None:
        ori = 'FF'

    genome = pysam.Fastafile("/Users/Nick_curie/Documents/Curie/Data/Genomes/Dmel_v6.12/Dmel_6.12.fasta")

    run_script(pos, n, split_read, genome, upstream_seq, downstream_seq, ori)

def reversed_seq(seq):
    return seq.translate(maketrans('ACGTacgt', 'TGCAtgca'))[::-1]


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
        position += 1
        i -= 1
    return(position , longest_hom, mhseq)


def getMechanism(homlen, inslen, delsize, templen):
    if inslen >= 10:
        mechanism="FoSTeS"
    elif templen >= 5:
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
    # print("Seq1:%s") % (seq1)
    # print("Seq2:%s") % (seq2)
    # print(seq1_start, seq1_end, seq2_start, seq2_end, seq)
    return(seq1_start, seq1_end, seq2_start, seq2_end, seq)

def run_script(pos, n, split_read, genome, upstream_seq, downstream_seq, ori):
    if not upstream_seq:
        (chrom1, bp1, bp2) = get_parts(pos)

        # bp1 -= 1
        bp2 -= 1
        upstream = bp1 - 200
        downstream = bp2 + 200

        upstream_seq = genome.fetch(chrom1, upstream, bp1)
        downstream_seq = genome.fetch(chrom1, bp2, downstream)
        # upstream_seq = reversed_seq(upstream_n)

    ##################
    ## Microhomology ##
    ##################

    print("Upstream:   %s") % (upstream_seq[-50:])
    print("Downstream: %s") % (downstream_seq[:50])

    if microhomology(upstream_seq, downstream_seq):
        (position, longest_hom, mhseq) =  microhomology(upstream_seq, downstream_seq)
    else:
        (position, longest_hom, mhseq) = (0,0,'')

    if longest_hom > 0:
        print("\n* Microhomology at breakpoint: %s (%s bp)") % (mhseq, longest_hom)
        downstream_spacer = (len(upstream_seq[-50:]) - len(mhseq))
        dmarker = " "*(downstream_spacer) + "^"*len(mhseq)
        downstream_buff = " "*downstream_spacer

        # upstream_base_end = bp1 - len(mhseq)
        # downstream_base_end = bp2 + len(mhseq)
        marker = "^"*len(mhseq)
        print(" Upstream:      %s") % (upstream_seq[-50:])
        print(" Downstream:    %s%s") % (downstream_buff, downstream_seq[:50])
        print(" Microhomology: %s\n") % (dmarker)

    else:
        print("\n* No microhomology found")
        longest_hom = 0

    if n is None:
        if longest_hom<=2:
            n=20
        else:
            n=50

    ##############
    ## Homology ##
    ##############

    if n >= len(upstream_seq):
         n = len(upstream_seq)

    # upstream_seq[-n:] takes n characters from the end of the string
    (downstream_start, downstream_end, upstream_start, upstream_end, homseq) = longestMatch(downstream_seq[0:n],upstream_seq[-n:])

    if len(homseq) > 0:
        # downstream_spacer = len(upstream_seq) - n  +  upstream_start - downstream_start
        print("* Longest homologous sequence +/- %s bp from breakpoint: %s (%s bp)\n") % (n, homseq, len(homseq))

        # What's this doing?
        seq_diff = len(upstream_seq[-n:]) - len(upstream_seq[-n:])

        upshift=0
        downshift=0

        if upstream_start - downstream_start > 0:
            downshift = upstream_start - downstream_start
        else:
            upshift = downstream_start - upstream_start

        upspace   = " "*(upshift+5)
        downspace = " "*downshift

        umarker = " "*(upstream_start+seq_diff) + "^"*len(homseq)
        dmarker = " "*(downstream_start+5) + "^"*len(homseq)

        # padded_downstream = " "*downstream_spacer + downstream_seq

        print(" Upstream:      %s%s--/--")   % (upspace, upstream_seq[-n:])
        print(" Homology:      %s%s")   % (upspace, umarker)
        print(" Downstream:    %s--/--%s")   % (downspace, downstream_seq[0:n])
        print(" Homology:      %s%s\n") % (downspace, dmarker)

    else:
        print("\n* No homology found +/- %s bp from breakpoint\n") % (n)


    deletion_size = 0
    deleted_bases = None
    insertion_size = 0
    inserted_seq = None
    templated_insertion_size = 0
    templated_insertion_seq = ''

    if split_read is not None:

        if ori == 'FF':
            print("Split read: %s\n") % (split_read)

        ###############
        ## Inserions ##
        ###############

        print("* Split read aligned to upstream sequence:")
        # Align split read to upstream seq (normal orientation)
        (upstream_start, upstream_end, split_start_up, split_end_up, upseq) = longestMatch(upstream_seq, split_read)

        ### Needs work. Currently, different output produced using consensus (longer than split). Should trom split reads if this is the case...
        if  upstream_start < split_start_up:
            print("here!!")

        ####

        aligned_up = len(upseq)
        deletion_size = len(upstream_seq[upstream_end:])
        deltype = 'up'

        if ori == 'FR':
            rev_split_read = reversed_seq(split_read[split_end_up:])
            print("Reversing downstream alignment of split read: %s => %s\n") % (split_read[split_end_up:], rev_split_read)
            (downstream_start, downstream_end, split_start, split_end, downseq) = longestMatch(downstream_seq, rev_split_read)

        else:
            (downstream_start, downstream_end, split_start, split_end, downseq) = longestMatch(downstream_seq, split_read[split_end_up:])

        # Calculate length of aligned sequences
        split_len = len(split_read)
        aligned_up = len(upseq)
        aligned_down = len(downseq)

        # Split read length - aligned portion = insetion size
        # insertion_size = int(split_len - aligned_up - aligned_down)

        insertion_size = split_start

        #######################
        #### Microhomology ####
        #######################

        if longest_hom > 0:

            if deletion_size == 0:
                deletion_size = downstream_start - longest_hom
                if deletion_size <= 0:
                    deletion_size = 0
                    deleted_bases = None
                if deletion_size > 0:
                    deltype = 'down'

            # microhomology with deletion
            if deletion_size > 0:
                if deltype == 'up':
                    deleted_bases = upstream_seq[-deletion_size:]

                    deletion_fill = "."*(deletion_size)
                    print("* Deletion of %s bases (%s)" % (deletion_size, upstream_seq[-deletion_size:]))

                    print(" Upstream:      %s") % (upstream_seq[upstream_start:])
                    print(" Split read:    %s%s--/--%s\n") % (split_read[split_start_up:split_end_up], deletion_fill, split_read[split_end_up:len(split_read)])

                    if ori == 'FF':
                        (downstream_start, downstream_end, split_start, split_end, downseq) = longestMatch(downstream_seq, split_read)
                        add2split = "*"*len(mhseq)
                        # print(split_read[split_end-add2split:split_end])
                        difference = (split_end_up - longest_hom) + deletion_size
                        seqbuffer = " "*(5+len(split_read[0:split_end_up-len(add2split)]))

                        print(" Split read:    %s--/--%s%s%s") % (split_read[0:split_end_up-len(add2split)], add2split, deletion_fill, split_read[split_end_up:])
                        print(" Downstream:    %s%s\n") % (seqbuffer, downstream_seq[0:downstream_end])
                    else:
                        (downstream_start, downstream_end, split_start, split_end, downseq) = longestMatch(downstream_seq, rev_split_read)

                        add2split = "*"*len(mhseq)
                        difference = (split_end_up - longest_hom) + deletion_size
                        seqbuffer = " "*(5+len(rev_split_read[0:split_end_up-len(add2split)]))

                        print(" Split read:    %s--/--%s%s%s") % (rev_split_read[0:split_end_up-len(add2split)], add2split, deletion_fill, rev_split_read[split_end_up:])
                        print(" Downstream:    %s%s\n") % (seqbuffer, downstream_seq[0:downstream_end])

                else:
                    deleted_bases = downstream_seq[:deletion_size]

                    deletion_fill = "."*(deletion_size)
                    print("* Deletion of %s bases (%s)" % (deletion_size, upstream_seq[-deletion_size:]))

                    print(" Upstream:      %s") % (upstream_seq[upstream_start:])
                    print(" Split read:    %s--/--%s\n") % (split_read[split_start_up:split_end_up], split_read[split_end_up:len(split_read)])

                    if ori == 'FF':
                        (downstream_start, downstream_end, split_start, split_end, downseq) = longestMatch(downstream_seq, split_read)
                    else:
                        (downstream_start, downstream_end, split_start, split_end, downseq) = longestMatch(downstream_seq, rev_split_read)
                    # if longest_hom > 0 and split_read[split_start_up:split_end_up] == mhseq[:-deleted_bases]:

                    # difference = downstream_start - split_start
                    # print(difference)

                    add2split = "*"*len(mhseq)
                    # print(split_read[split_end-add2split:split_end])
                    difference = (split_end_up - longest_hom) + deletion_size
                    seqbuffer = " "*(5+len(split_read[0:split_end_up-len(add2split)]))

                    print(" Split read:    %s--/--%s%s%s") % (split_read[0:split_end_up-len(add2split)], deletion_fill, add2split, split_read[split_end_up:])
                    print(" Downstream:    %s%s\n") % (seqbuffer, downstream_seq[0:downstream_end])

            # microhomology with no deletion
            else:
                deleted_bases = None
                print("\n* No deletion at breakpoint")

                print(" Upstream:      %s") % (upstream_seq[upstream_start:])
                print(" Split read:    %s--/--%s\n") % (split_read[split_start_up:split_end_up], split_read[split_end_up:len(split_read)])

                # microhomology with insetion
                if insertion_size > 0:

                    if ori == 'FF':
                        (downstream_start, downstream_end, split_start, split_end, downseq) = longestMatch(downstream_seq, split_read[split_end_up:])
                        inserted_seq = split_read[split_end_up:split_end_up+insertion_size]
                        seqbuffer = " "*((aligned_up+5+insertion_size+2+2))
                        add2split = "*"*len(mhseq)
                        add2split = ''

                        print(" Split read:    %s--/--[%s]--%s%s") % (split_read[0:split_end_up], split_read[split_end_up:split_end_up+insertion_size],add2split, split_read[split_end_up+insertion_size:])
                        print(" Downstream:    %s%s\n") % (seqbuffer, downstream_seq)
                    else:
                        (downstream_start, downstream_end, split_start, split_end, downseq) = longestMatch(downstream_seq, rev_split_read)
                        inserted_seq = rev_split_read[split_start-insertion_size:split_start]

                        print("* %s bp insertion '%s' at breakpoint\n") % (insertion_size, inserted_seq)

                        seqbuffer = " "*((aligned_up+5+insertion_size+2+2))
                        add2split = "*"*len(mhseq)
                        add2split = ''

                        print(" Split read:    %s--/--[%s]--%s%s") % (rev_split_read[0:split_end_up], inserted_seq, add2split, rev_split_read[split_start:split_end])
                        print(" Downstream:    %s%s\n") % (seqbuffer, downstream_seq[:split_end-insertion_size])
                # microhomology with no insetion and no deletion
                else:

                    if ori == 'FF':
                        (downstream_start, downstream_end, split_start, split_end, downseq) = longestMatch(downstream_seq, split_read[split_end_up:])
                        difference = (downstream_start-split_start)

                        seqbuffer = " "*(5+len(split_read[0:split_end_up]))
                        nonaligned = "*"*difference

                        print(" Split read:    %s--/--%s%s") % (split_read[0:split_end_up], nonaligned, split_read[split_end_up:])
                        print(" Downstream:    %s%s\n") % (seqbuffer, downstream_seq[0:downstream_end])
                    else:
                        (downstream_start, downstream_end, split_start, split_end, downseq) = longestMatch(downstream_seq, rev_split_read[split_start:])
                        difference = (downstream_start-split_start)

                        seqbuffer = " "*(5+len(split_read[0:split_end_up]))
                        nonaligned = "*"*difference

                        print(" Split read:    %s--/--%s%s") % (split_read[0:split_end_up], nonaligned, rev_split_read[split_start:])
                        print(" Downstream:    %s%s\n") % (seqbuffer, downstream_seq[0:downstream_end])


        ##########################
        #### No Microhomology ####
        ##########################

        else:

            if deletion_size == 0:
                deletion_size = downstream_start
                deltype = 'down'

            # Deletion without microhomology
            if deletion_size > 0:
                if deltype == 'up':
                    deleted_bases = upstream_seq[-deletion_size:]
                else:
                    deleted_bases = downstream_seq[:deletion_size]

                # deleted_bases = upstream_seq[-deletion_size:]
                deletion_fill = "."*(deletion_size)
                print("* Deletion of %s bases (%s)" % (deletion_size, deleted_bases))

                if deltype == 'up':
                    print(" Upstream:      %s") % (upstream_seq[upstream_start:])
                    print(" Split read:    %s%s--/--%s\n") % (split_read[split_start_up:split_end_up], deletion_fill, split_read[split_end_up:len(split_read)])

                    (downstream_start, downstream_end, split_start, split_end, downseq) = longestMatch(downstream_seq, split_read)

                    difference = (split_start - downstream_start)
                    seqbuffer = " "*(5+len(split_read[0:split_end_up]))

                    print(" Split read:    %s--/--%s") % (split_read[0:split_end_up], split_read[split_start:])
                    print(" Downstream:    %s%s\n") % (seqbuffer, downstream_seq[0:downstream_end])
                else:
                    print(" Upstream:      %s") % (upstream_seq[upstream_start:])
                    print(" Split read:    %s--/--%s\n") % (split_read[split_start_up:split_end_up], split_read[split_end_up:len(split_read)])

                    (downstream_start, downstream_end, split_start, split_end, downseq) = longestMatch(downstream_seq, split_read)

                    difference = (split_start - downstream_start)
                    seqbuffer = " "*(5+len(split_read[0:split_end_up]))

                    print(" Split read:    %s--/--%s%s") % (split_read[0:split_end_up], deletion_fill,split_read[split_start:])
                    print(" Downstream:    %s%s\n") % (seqbuffer, downstream_seq[0:downstream_end])

            # Insertion without microhomology
            if insertion_size > 0:
                if ori == 'FF':
                    (downstream_start, downstream_end, split_start, split_end, downseq) = longestMatch(downstream_seq, split_read[split_end_up:])
                else:
                    (downstream_start, downstream_end, split_start, split_end, downseq) = longestMatch(downstream_seq, rev_split_read)
                inserted_seq = split_read[split_end_up:split_end_up+insertion_size]
                print("* %s bp insertion '%s' at breakpoint\n") % (insertion_size, inserted_seq)

                # This is lifted from below - not fully tested
                print(" Upstream:      %s") % (upstream_seq[upstream_start:])
                print(" Split read:    %s--/--%s\n") % (split_read[split_start_up:split_end_up], split_read[split_end_up:])

                seqbuffer = " "*((aligned_up+5+insertion_size+2+2))

                # What's this doing?
                add2split = "*"*len(mhseq)
                # add2split = ''

                print(" Split read:    %s--/--[%s]--%s%s") % (split_read[:split_end_up], split_read[split_end_up:split_end_up+insertion_size],add2split, split_read[split_end_up+insertion_size:])
                print(" Downstream:    %s%s\n") % (seqbuffer, downstream_seq[:split_end])

            if insertion_size == 0 and deletion_size == 0:
                print("\n* No microhomolgy found, and no deletion or insertion at breakpoint")

                print(" Upstream:      %s") % (upstream_seq[upstream_start:])
                print(" Split read:    %s--/--%s\n") % (split_read[split_start_up:split_end_up], split_read[split_end_up:len(split_read)])

                if ori == 'FF':
                    (downstream_start, downstream_end, split_start, split_end, downseq) = longestMatch(downstream_seq, split_read)

                    seqbuffer = " "*(5+len(split_read[0:split_end_up]))
                    # nonaligned = "-"*difference
                    print(" Split read:    %s--/--%s") % (split_read[0:split_end_up], split_read[split_end_up:])
                    print(" Downstream:    %s%s\n") % (seqbuffer, downstream_seq[0:split_end])

                else:
                    (downstream_start, downstream_end, split_start, split_end, downseq) = longestMatch(downstream_seq, rev_split_read)

                    seqbuffer = " "*(5+len(split_read[0:split_end_up]))
                    # nonaligned = "-"*difference
                    print(" Split read:    %s--/--%s") % (split_read[0:split_end_up], rev_split_read[split_start:])
                    print(" Downstream:    %s%s\n") % (seqbuffer, downstream_seq[0:split_end])


        ########################
        #  Realign insertion ###
        ########################

        templated_up = ''
        templated_down = ''
        if insertion_size >= 3:
            # n=20
            (inserted_start, inserted_end, upstream_start, upstream_end, aligned) = longestMatch(inserted_seq,upstream_seq[-n:])
            if len(aligned) >= 3:
                templated_up = aligned
                templated_insertion_seq = templated_up
                templated_insertion_size = len(templated_up)
                splitbuffer = " "*upstream_start
                insertion_pos = (len(upstream_seq[-n:]) - upstream_end)
                print("* %s bp of inserted sequence -%s bps from breakpoint on upstream sequence") % (len(templated_up),insertion_pos)

                print(" Upstream:      %s--/--") % (upstream_seq[-n:])
                print(" Insertion:     %s%s\n") % (splitbuffer, templated_up)
            else:
                print("Could not find at least 3 bases of inserted sequence in upstream region")
            (downstream_start, downstream_end, inserted_start, inserted_end, aligned) = longestMatch(downstream_seq[:n], inserted_seq)

            if len(aligned) >= 3:
                templated_down = aligned
                templated_insertion_seq = templated_down
                if len(templated_down) > templated_insertion_size:
                    templated_insertion_size = len(templated_down)

                splitbuffer = " "*(downstream_start+5)
                print("* %s bp of inserted sequence +%s bps from breakpoint on downstream sequence") % (len(templated_down),downstream_start)
                print(" Downstream:     --/--%s") % (downstream_seq[:n])
                print(" Insertion:      %s%s\n") % (splitbuffer, templated_down)
            else:
                print("Could not find at least 3 bases of inserted sequence in downstream region")


         # else:
         #      print("No split read sequence provided. Unable to find insetion at bp")
         #      insertion_size = 0
    #
    # ###############
    # ## Mechanism ##
    # ###############

    # Calculate mechanism
    mechanism=getMechanism(longest_hom, insertion_size, deletion_size, templated_insertion_size)
    print("* Mechanism: %s") % (mechanism)
    print("  * %s bp '%s' microhomology at breakpoints") % (longest_hom, mhseq)
    print("  * %s bp deletion") % (deletion_size)
    print("  * %s bp insertion") % (insertion_size)
    if insertion_size > 3:
        print("  * %s bp of inserted seq '%s' found locally templated") % (templated_insertion_size, templated_insertion_seq)


    return(longest_hom, mhseq, homseq, deletion_size, deleted_bases, insertion_size, inserted_seq, templated_up, templated_down, templated_insertion_size)



    #



if __name__ == '__main__':
    main()
