from Bio.Align.Applications import ClustalOmegaCommandline


upstream_seq = 'AAAAAAAAAACCTTG'
downstream_seq = 'CTTGGCATTACGTA'
split_read = 'AACCCTTGGCATTACGT'

tempfile_up   = 'data/up_split.fa'
tempfile_down = 'data/down_split.fa'

with open(tempfile_up, 'w') as up_tmp, open(tempfile_down, 'w') as down_tmp:
    up_tmp.write(">upstream\n%s\n" % upstream_seq)
    up_tmp.write(">split\n%s\n" % split_read)

    down_tmp.write(">downstream\n%s\n" % downstream_seq)
    down_tmp.write(">split\n%s\n" % split_read)



temp_aligned_up = 'data/temp_aligned_up'
clustalomega_cline = ClustalOmegaCommandline(infile=tempfile_up, outfile=temp_aligned_up, verbose=True, auto=True,force=True, outfmt="st")

print(clustalomega_cline)

clustalomega_cline()

from Bio import AlignIO
alignment = AlignIO.read(open(temp_aligned_up), "stockholm")
print("Alignment length %i" % alignment.get_alignment_length())

for record in alignment :
    print(record.seq + " " + record.id)
