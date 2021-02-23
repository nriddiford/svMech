#!/usr/bin/env python
import os, re, sys

from collections import defaultdict
from optparse import OptionParser

import pysam
import pandas as pd


def parse_variants(options):
    df = pd.read_csv(options.in_file, delimiter = "\t")
    df = df.where((pd.notnull(df)), None)

    df.columns = ["event", "bp_no", "sample", "genotype", "chrom", "bp", "gene", "feature", "chrom2", "bp2", "gene2", "feature2", "type", "length"]
    filtered_df = pd.DataFrame()

    seen_events = defaultdict(int)

    exclude_sample = ["A373R1", "A373R7", "A512R17", "A373R11", "A785-A788R1", "A785-A788R11", "A785-A788R3", "A785-A788R5", "A785-A788R7", "A785-A788R9"]

    filtered_df = df[~df['sample'].isin(exclude_sample)]
    filtered_df = filtered_df[filtered_df['genotype'] == 'somatic_tumour']
    filtered_df = filtered_df[filtered_df['type'] == 'DEL']


    for index, variant in filtered_df.iterrows():
        key = '_'.join([variant['sample'], str(variant['event']), str(variant['bp_no']), str(variant['length'])])
        seen_events[key] += 1

        if seen_events[key] > 1:
            print("Seen this event before: %s, %s") % (variant['sample'], str(variant['event']))
            print seen_events[key]
            continue

    return filtered_df


def parseLoci(options):
    df = pd.read_csv(options.in_file, delimiter="\t")
    df.columns = ["chrom", "start", "end"]

    df.rename(columns={'start': 'bp'}, inplace=True)

    df['sample'] = 'sim'
    # df['type'] = df.groupby(['sample']).ngroup()
    df['type'] = df.reset_index().index

    return df


def parseBed(options):
    df = pd.read_csv(options.in_file, delimiter="\t")
    df.columns = ["chrom", "start", "end"]

    df['sample'] = 'sim'
    # df['type'] = df.groupby(['sample']).ngroup()
    df['type'] = df.reset_index().index

    return df


def parse_SNV(options):
    df = pd.read_csv(options.in_file, delimiter="\t")
    df.columns = ["sample", "chrom", "pos", "ref", "alt", "tri", "trans", "decomposed_tri", "grouped_trans", "a_freq", "caller", "feature", "gene", "id"]

    df.rename(columns={'pos': 'bp', 'grouped_trans': 'type'}, inplace=True)

    filtered_df = pd.DataFrame()

    exclude_sample = ["A373R1", "A373R7", "A512R17", "A373R11", "A785-A788R1", "A785-A788R11", "A785-A788R3", "A785-A788R5", "A785-A788R7", "A785-A788R9"]

    filtered_df = df[~df['sample'].isin(exclude_sample)]
    filtered_df = filtered_df[filtered_df['a_freq']>=0.3]

    return filtered_df


def getSeqs(df, options):

    window = int(options.distance)

    genome = pysam.Fastafile(options.genome)

    with open(options.out_file, 'w+') as seqs_out:
        for index, variant in df.iterrows():

            upstream = variant['bp'] - window/2
            downstream = variant['bp'] + window/2

            surrounding_seq = genome.fetch(variant['chrom'], upstream, downstream)

            # downstream_seq = genome.fetch(variant['chrom'], variant['bp'], downstream)

            if options.type == 'bp':
                seqs_out.write(">%s_%s_%s_%s_%s\n%s\n" % (variant['sample'], variant['type'], variant['length'], variant['bp_no'], window, surrounding_seq))
                # seqs_out.write(">%s_%s_%s_%s_%s_%s\n%s\n" % (variant['sample'], variant['type'], variant['length'], variant['bp_no'], "Downstream", window, downstream_seq))
            elif options.type == 'snv':
                # seqs_out.write(">%s_%s_%s:%s_%s_%s_%s\n%s\n" % (variant['sample'], variant['type'], variant['chrom'], variant['bp'], variant['a_freq'], "Upstream", window, upstream_seq))
                seqs_out.write(">%s_%s_%s:%s_%s_%s\n%s\n" % (variant['sample'], variant['type'], variant['chrom'], variant['bp'], variant['a_freq'], window, surrounding_seq))
            elif options.type == 'locus':
                seqs_out.write(">%s_%s_%s:%s_%s\n%s\n" % (variant['sample'], variant['type'], variant['chrom'], variant['bp'], window, surrounding_seq))



def get_args():
    parser = OptionParser()

    parser.add_option("-i",
                      "--in_file",
                      dest = "in_file",
                      action = "store",
                      help = "A tsv file containing variant information ",
                      metavar = "FILE")

    parser.add_option("-t",
                      "--type",
                      dest = "type",
                      action = "store",
                      help = "Type of variant file [bp, snv, locus or bed]")

    parser.add_option("-g",
                      "--genome",
                      dest="genome",
                      action="store",
                      help="Genome fasta file",
                      metavar="FILE")


    parser.add_option("-o",
                      "--out_file",
                      dest = "out_file",
                      action = "store",
                      help = "File to output breakpoint fasta",
                      metavar ="FILE")

    parser.add_option("-d",
                      "--distance",
                      dest = "distance",
                      action = "store",
                      help = "Distance from breakpoint to return")


    parser.set_defaults(in_file = 'data/all_somatic_tumour.tsv',
                        distance = 100,
                        type = 'bp',
                        genome = "/Users/Nick_curie/Documents/Curie/Data/Genomes/Dmel_v6.12/Dmel_6.12.fasta")

    options, args = parser.parse_args()

    if not options.out_file:
        dist = str(options.distance)
        filename = 'somatic_breakpoints' + '_' + dist + '.fa'
        options.out_file = os.path.join('out', filename)



    if options.in_file is None:
        parser.print_help()
        print
    return (options, args)


def main():
    options, args = get_args()

    if not os.path.exists(os.path.dirname(options.out_file)):
        try:
            os.makedirs(os.path.dirname(options.out_file))
        except OSError as exc:  # Guard against race condition
            if exc.errno != errno.EEXIST:
                raise

    if options.in_file and options.distance:
        try:
            if options.type == 'bp':
                print("Parsing breakpoints")
                vars = parse_variants(options)
            elif options.type == 'snv':
                print("Parsing SNVs")
                vars = parse_SNV(options)
            elif options.type == 'bed':
                print("Parsing bed file")
                vars = parseBed(options)
            elif options.type == 'locus':
                print("Parsing genomic loci")
                vars = parseLoci(options)

            getSeqs(vars, options)

        except IOError as err:
            sys.stderr.write("IOError " + str(err) + "\n")
            return


if __name__ == "__main__":
    sys.exit(main())
