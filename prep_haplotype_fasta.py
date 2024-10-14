from Bio import SeqIO
import sys
import re

# Process output of assembly-to-fasta to prepare it for hic file generation. 
# Takes care of:
# - Renaming fasta sequences so that they are unique to each haplotype and 
#   differentiate chromosomes (LG, linkage groups) from unincorporated scaffolds
# - Removing trailing NNNNN from sequences

# Example use:
# python ../prep_haplotype_fasta.py NBI_hap1_v0.5.fasta NBI_hap1_v0.5.renamed.fasta NBIv0.5_hap1 8
# Where 8 signififies that the first 8 scaffolds are linkage groups (chromosomes)

in_name = sys.argv[1]
out_name = sys.argv[2]
prefix = sys.argv[3]
lg = int(sys.argv[4])

renumber = True
if len(sys.argv)>5:
  if sys.argv[5] == 'r':
    print("Renumbering scaffolds consecutively")
    renumber=True

# Example input names:
# >Scaffold1
# >Scaffold2
# >Scaffold3

# Example output names:
# >NBIv0.4_hap1_LG8
# >NBIv0.4_hap1_scaf76

def rename( header ):
  number = int(header[8:])
  if number <= lg:
    return( f"{prefix}_LG{number}" )
  else:
    return( f"{prefix}_scaf{number}" )

def rename_renumber( number ):
  if number <= lg:
    return( f"{prefix}_LG{number}" )
  else:
    return( f"{prefix}_scaf{number}" )


print("id length N_start N_end")

with open(out_name, "w") as output_handle:
  i = 1
  for record in SeqIO.parse(in_name, "fasta"):
    if renumber:
      record.id = rename_renumber( i )
    else:
      record.id = rename( record.id )
    record.description = ""
    
    N_start = False
    N_end = False
    
    if record.seq[0] == "N":
      N_start = True
      record.seq = re.sub("^N+", "", record.seq)
    
    if record.seq[-1] == "N":
      N_start = True
      record.seq = re.sub("N+$", "", record.seq)
    
    print(f"{record.id} {len(record.seq)} {N_start} {N_end}")
    SeqIO.write(record, output_handle, "fasta")
    i = i + 1


