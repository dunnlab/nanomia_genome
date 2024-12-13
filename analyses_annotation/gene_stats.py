import argparse
import pandas as pd


def main(input_file, prefix):
    df_gff = pd.read_csv(input_file, sep="\t", comment="#", header=None)
    df_gff.columns = ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]

    df_gff = df_gff.sort_values(by=["seqid", "start"])

    if prefix:
        df_gff = df_gff[df_gff["seqid"].str.startswith(prefix)].copy()

    print(f"Number of rows read: {len(df_gff)}")
    if prefix:
        print(f"Number of rows after filtering by prefix: {len(df_gff)}")

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

    df_gff['length'] = df_gff['end'] - df_gff['start'] + 1

    # Create gene data frame
    df_genes = df_gff[df_gff["type"] == "gene"].copy()
    df_genes["mRNA_ID"] = df_genes["ID"] + "-RA"

    exon_counts = df_gff[df_gff["type"] == "exon"].groupby("Parent").size()
    df_genes["num_exons"] = df_genes["mRNA_ID"].map(exon_counts).fillna(0).astype(int)

    exon_lengths = df_gff[df_gff["type"] == "exon"].groupby("Parent")["length"].sum()
    df_genes["exons_length"] = df_genes["mRNA_ID"].map(exon_lengths).fillna(0).astype(int)
    df_genes["introns_length"] = df_genes["length"] - df_genes["exons_length"]#

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
    df_introns["intron_length"] = df_introns.groupby("Parent")["end"].shift(1).fillna(0)
    df_introns["intron_length"] = df_introns["start"] - df_introns["intron_length"] - 1
    df_introns["intron_length"] = df_introns["intron_length"].astype(int)

    # Keep only rows where the Parent is the same as the previous row
    df_introns = df_introns[df_introns["Parent"] == df_introns["Parent"].shift(1)].copy()
    df_introns["type"] = "intron"

    df_introns = df_introns[["seqid", "type", "intron_length", "ID", "Parent"]]
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
