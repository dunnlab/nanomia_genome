from Bio import SeqIO
from Bio.Seq import Seq
import sys

forward = "TTAGGG"
reverse = str(Seq(forward).reverse_complement())

# Look for repeats, it is what we are interested in and there are many fewer than
# there are singletons
forward2 = forward + forward
reverse2 = reverse + reverse

in_name = sys.argv[1]

print("scaffold sequence position")

for record in SeqIO.parse(in_name, "fasta"):
  name = record.id
  seq = record.seq
  target_length = len(forward2)
  for i in range(0, len(seq)-target_length-1):
    kmer = seq[i: (i + target_length) ]
    if (kmer == forward2):
      print(name, forward, i)
    elif (kmer == reverse2):
      print(name, reverse, i)