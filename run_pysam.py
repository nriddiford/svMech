import pysam
samfile = pysam.AlignmentFile("data/R5_Del.bam", "rb")

iter = samfile.fetch("X", 8812326,8815425)
for x in iter:
    print (str(x))
