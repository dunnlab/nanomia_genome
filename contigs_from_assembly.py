# From https://github.com/conchoecia/genome_assembly_pipelines/blob/f9eba44569b0f1c8686268f24650728006d118e4/scripts/GAP_hic_map7#L659
import sys
from Bio import SeqIO


assem = sys.argv[1]

outbed = sys.argv[2]
# this block of code from https://www.biostars.org/p/133742/
#import the SeqIO module from Biopython
outhandle = open(outbed, "w")
with open(assem, mode="r") as fasta_handle:
    for record in SeqIO.parse(fasta_handle, "fasta"):
        start_pos=0
        counter=0
        gap=False
        gap_length = 0
        for char in record.seq:
            if char == 'N':
                if gap_length == 0:
                    gap_length = 1
                    gap = True
                    print("{} {} {}".format(
                        record.id,
                        start_pos,
                        counter),
                          file = outhandle)
                else:
                    gap_length += 1
            else:
                if gap:
                    gap_length = 0
                    gap = False
                    start_pos=counter
            counter += 1
        print("{} {} {}".format(
            record.id,
            start_pos,
            counter),
              file = outhandle)
outhandle.close()