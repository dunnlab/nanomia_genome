from Bio import SeqIO
import sys
import re


in_name = sys.argv[1]

for record in SeqIO.parse(in_name, "fasta"):
  name = record.id
  seq = record.seq

  # Count the number of poly Ns
  n = len(re.findall("N+", str(seq))) + 1

  # Get the length of seq without the poly Ns
  len_nogaps = len(re.sub("N+", "", str(seq)))

  print( name, len(seq),len_nogaps,  n, sep="\t" )