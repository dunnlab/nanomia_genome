import argparse
import pandas as pd
from enum import Enum


def gene_stats(input_file, prefix):

    class FileType(Enum):
        GFF = "GFF"
        GTF = "GTF"

    if input_file.endswith(".gff"):
        file_type = FileType.GFF
    elif input_file.endswith(".gtf"):
        file_type = FileType.GTF
    else:
        raise ValueError("File type not supported")

    print(f"Reading file: {input_file}")
    print(f"File type: {file_type}")

    df_gff = pd.read_csv(input_file, sep="\t", comment="#", header=None)
    df_gff.columns = ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]

    df_gff = df_gff.sort_values(by=["seqid", "start"])

    if prefix:
        df_gff = df_gff[df_gff["seqid"].str.startswith(prefix)].copy()

    print(f"Number of rows read: {len(df_gff)}")

    df_gff.drop_duplicates(inplace=True)
    print(f"Number of rows after removing duplicates: {len(df_gff)}")

    if prefix:
        print(f"Number of rows after filtering by prefix: {len(df_gff)}")

    df_gff['gene_ID'] = None

    # Handle the gff case
    if file_type == FileType.GFF:
        # Extract attributes in a vectorized way
        df_gff['ID'] = df_gff['attributes'].str.extract(r'ID=([^;]+)')
        df_gff['Parent'] = df_gff['attributes'].str.extract(r'Parent=([^;]+)')
        df_gff['Note'] = df_gff['attributes'].str.extract(r'Note=([^;]+)')
        
        # Assign gene_ID vectorized
        df_gff['gene_ID'] = None
        df_gff.loc[df_gff["type"] == "gene", "gene_ID"] = df_gff["ID"]
        df_gff.loc[df_gff["type"] == "mRNA", "gene_ID"] = df_gff["Parent"]
        df_gff.loc[df_gff["type"] == "CDS", "gene_ID"] = df_gff["Parent"]
        df_gff.loc[df_gff["type"] == "exon", "gene_ID"] = df_gff["Parent"].str.replace("-RA", "", regex=False)

    # Handle the gtf case
    elif file_type == FileType.GTF:
        # GTF-specific attribute parsing
        df_gff['gene_ID'] = df_gff['attributes'].str.extract(r'gene_id "([^"]+)"')
        df_gff['transcript_ID'] = df_gff['attributes'].str.extract(r'transcript_id "([^"]+)"')
        
        # Ensure gene_ID is not empty for rows where it's required
        df_gff.loc[df_gff['type'] == 'gene', 'gene_ID'] = df_gff['gene_ID']

    else:
        raise ValueError("File type not supported")

    df_gff['length'] = df_gff['end'] - df_gff['start'] + 1

    # Create gene data frame
    df_genes = df_gff[df_gff["type"] == "gene"].copy()

    exon_counts = df_gff[df_gff["type"] == "exon"].groupby("gene_ID").size()
    df_genes["num_exons"] = df_genes["gene_ID"].map(exon_counts).fillna(0).astype(int)

    exon_lengths = df_gff[df_gff["type"] == "exon"].groupby("gene_ID")["length"].sum()
    df_genes["exons_length"] = df_genes["gene_ID"].map(exon_lengths).fillna(0).astype(int)
    df_genes["introns_length"] = df_genes["length"] - df_genes["exons_length"]

    df_genes["intron_ratio"] = df_genes["introns_length"] / df_genes["exons_length"]

    # Add column intergenic_distance that is the difference between the start of the current
    # gene and the end of the previous gene
    df_genes["intergenic_distance"] = df_genes["start"].shift(1).fillna(0)
    df_genes["intergenic_distance"] = df_genes["start"] - df_genes["intergenic_distance"] - 1

    # Set intergenic_distance to None for the first gene in each sequence
    df_genes.loc[df_genes["seqid"] != df_genes["seqid"].shift(1), "intergenic_distance"] = None

    # Add column intergenic_length that is the difference between the start of the current
    # gene and the end of the previous gene
    df_genes["intergenic_length"] = df_genes["start"].shift(1).fillna(0)
    df_genes["intergenic_length"] = df_genes["start"] - df_genes["intergenic_length"] - 1

    # Create exon data frame
    df_exons = df_gff[df_gff["type"] == "exon"].copy()

    # Create intron data frame
    df_introns = df_exons.copy()

    # Add column intron_length that is the difference between the start of the current
    # exon and the end of the previous exon
    df_introns["intron_length"] = df_introns.groupby("gene_ID")["end"].shift(1).fillna(0)
    df_introns["intron_length"] = df_introns["start"] - df_introns["intron_length"] - 1
    df_introns["intron_length"] = df_introns["intron_length"].astype(int)

    # Keep only rows where the gene_ID is the same as the previous row
    df_introns = df_introns[df_introns["gene_ID"] == df_introns["gene_ID"].shift(1)].copy()
    df_introns["type"] = "intron"

    df_introns = df_introns[["seqid", "type", "intron_length", "gene_ID"]]
    df_introns.rename(columns={"intron_length": "length"}, inplace=True)

    return df_genes, df_exons, df_introns


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate gene statistics")
    parser.add_argument("input", help="Input gff file")
    parser.add_argument("-p", "--prefix", help="Prefix to filter seqid column")
    args = parser.parse_args()
    (df_genes, df_exons, df_introns) = gene_stats(args.input, args.prefix)

    print("Gene statistics")
    print(f" Number: {len(df_genes)}")
    print(f" Min length: {df_genes['length'].min()}")
    print(f" Max length: {df_genes['length'].max()}")
    print(f" Mean length: {df_genes['length'].mean()}")
    print(f" Median length: {df_genes['length'].median()}")
    print(f" Total length: {df_genes['length'].sum()}")

    print("Inter-genic statistics")
    intergenic_distance = df_genes["intergenic_distance"].dropna()
    print(f" Number: {len(intergenic_distance)}")
    print(f" Min length: {intergenic_distance.min()}")
    print(f" Max length: {intergenic_distance.max()}")
    print(f" Mean length: {intergenic_distance.mean()}")
    print(f" Median length: {intergenic_distance.median()}")
    print(f" Total length: {intergenic_distance.sum()}")

    print("Exon statistics")
    print(f" Number: {len(df_exons)}")
    print(f" Min length: {df_exons['length'].min()}")
    print(f" Max length: {df_exons['length'].max()}")
    print(f" Mean length: {df_exons['length'].mean()}")
    print(f" Median length: {df_exons['length'].median()}")
    print(f" Total length: {df_exons['length'].sum()}")

    print("Intron statistics")
    print(f" Number: {len(df_introns)}")
    print(f" Min length: {df_introns['length'].min()}")
    print(f" Max length: {df_introns['length'].max()}")
    print(f" Mean length: {df_introns['length'].mean()}")
    print(f" Median length: {df_introns['length'].median()}")
    print(f" Total length: {df_introns['length'].sum()}")

    print(f"Intron ratio: {df_introns['length'].sum() / df_exons['length'].sum():.3f}")
