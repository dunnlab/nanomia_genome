import sys
import pandas as pd

in_file = sys.argv[1]

# Parse GFF file
df = pd.read_csv(in_file, sep='\t', header=None, names=['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes'])

# Only retain rows with seqid that contains 'LG'
df = df[df['seqid'].str.contains('LG')]

# Get a list of unique seqids
seqids = df['seqid'].unique()

# Get a list of unique attributes
attributes = df['attributes'].unique()

for each seqid in seqids:
  for each attribute in attributes:
    # Get a list of all rows that have the seqid and attribute
		rows = df[(df['seqid'] == seqid) & (df['attributes'] == attribute)]
		# If there are more than one row, print the seqid and attribute
		if len(rows) > 1:
			print(seqid, attribute)




