import argparse
import pandas as pd
from enum import Enum


def main(input_file, prefix):

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

        # Extract attributes and handle missing values with an explicit loop
        attributes_dict = {'ID': [], 'Parent': [], 'Note': []}
        for attributes in df_gff['attributes']:
            attr_dict = {}
            for key_value in attributes.split(";"):
                if "=" in key_value:
                    key, value = key_value.split("=")
                    attr_dict[key] = value
            for key in attributes_dict:
                attributes_dict[key].append(attr_dict.get(key))

        df_gff['ID'] = attributes_dict['ID']
        df_gff['Parent'] = attributes_dict['Parent']
        df_gff['Note'] = attributes_dict['Note']

        # Loop over the rows and fill the gene_ID column
        # if gene, gene_ID = ID
        # if mRNA, gene_ID = Parent
        # if CDS, gene_ID = Parent
        # if exon, gene_ID = Parent with terminal -RA removed
        gene_ID = None
        for i, row in df_gff.iterrows():
            if row["type"] == "gene":
                gene_ID = row["ID"]
            elif row["type"] == "mRNA":
                gene_ID = row["Parent"]
            elif row["type"] == "CDS":
                gene_ID = row["Parent"]
            elif row["type"] == "exon":
                gene_ID = row["Parent"].replace("-RA", "")
            df_gff.loc[i, "gene_ID"] = gene_ID

    # Handle the gtf case
    elif file_type == FileType.GTF:
        # Extract attributes specific to GTF
        attributes_dict = {'gene_id': [], 'transcript_id': []}
        for attributes in df_gff['attributes']:
            attr_dict = {'gene_id': None, 'transcript_id': None}
            for entry in attributes.split(";"):
                entry = entry.strip()
                if entry:
                    key_value = entry.split(" ", 1)  # Split on the first space only
                    if len(key_value) == 2:
                        key, value = key_value
                        value = value.replace('"', '')
                        if key in attr_dict:
                            attr_dict[key] = value
            for key in attributes_dict:
                attributes_dict[key].append(attr_dict.get(key))

        # Create columns for GTF file
        df_gff['gene_ID'] = attributes_dict['gene_id']
        df_gff['transcript_ID'] = attributes_dict['transcript_id']

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
    (df_genes, df_exons, df_introns) = main(args.input, args.prefix)

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

    # print any lines from df_introns that have negative length
    print(df_introns[df_introns["length"] < 0])
