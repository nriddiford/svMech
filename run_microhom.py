from microhom import run_script
import pysam

deletions_file='data/deletion_input.txt'

# genome = pysam.Fastafile("/Users/Nick_curie/Documents/Curie/Data/Genomes/Dmel_v6.12/Dmel_6.12.fasta")
genome = pysam.Fastafile("/Users/Nick/Documents/Curie/Data/Genomes/Dmel_v6.12/Dmel_6.12.fasta")


# def run_script(pos, n, split_read, genome, upstream_seq, downstream_seq, ori, cut):
dels_out = open('data/deletions.txt', 'w')
headers = ["Sample", "Length", "microhomolgy", "extended_homology", "Deleted_bases", "Inserted_bases", "templated_up", "templated_down", "Mechanism"]
dels_out.write('\t'.join(headers) + '\n')

with open(deletions_file, 'r') as dels:
    for l in dels:
        parts = l.rstrip().split('\t')
        ori = 'FF'

        try:
            (sample, length, split_read, position, cut) = parts
            print(position, split_read, cut)

            (longest_hom, mhseq, homseq, deletion_size, deleted_bases, insertion_size, inserted_seq, templated_up, templated_down, templated_insertion_size, mechanism) = run_script(position, None, split_read, genome, '', '', ori, cut)
        except ValueError:
            cut = 200
            (sample, length, split_read, position) = parts
            print(position, split_read, cut)
            (longest_hom, mhseq, homseq, deletion_size, deleted_bases, insertion_size, inserted_seq, templated_up, templated_down, templated_insertion_size, mechanism) = run_script(position, None, split_read, genome, '', '', ori, 200)

        print("Sample: %s\nTemplated up: %s\nTemplated down: %s\n") % (sample, templated_up, templated_down)

        # outline = [sample, length, mhseq, homseq, deleted_bases, inserted_seq, templated_up, templated_down, "test"]
        # d.append((sample, length, mhseq, homseq, deleted_bases, inserted_seq, templated_up, templated_down, "test"))

        # dels_out.write('\t'.join(outline) + '\n')
        dels_out.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (sample, length, mhseq, homseq, deleted_bases, inserted_seq, templated_up, templated_down, mechanism))
# pd.DataFrame(d, columns=("Sample", "Length", "microhomolgy", "extended_homology", "Deleted_bases", "Inserted_bases", "templated_up", "templated_down", "Mechanism"))
