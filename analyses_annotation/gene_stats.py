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

    gene_ids = set(df_gff[df_gff["type"] == "gene"]["ID"])
    print(f"Number of genes: {len(gene_ids)}")

    df_genes = df_gff[df_gff["type"] == "gene"].copy()
    df_genes["mRNA_ID"] = df_genes["ID"] + "-RA"

    exon_counts = df_gff[df_gff["type"] == "exon"].groupby("Parent").size()
    df_genes["num_exons"] = df_genes["mRNA_ID"].map(exon_counts).fillna(0).astype(int)

    # Create a new column with the sum of the lengths of all exons for each gene
    exon_lengths = df_gff[df_gff["type"] == "exon"].groupby("Parent")["length"].sum()
    df_genes["exons_length"] = df_genes["mRNA_ID"].map(exon_lengths).fillna(0).astype(int)
    df_genes["introns_length"] = df_genes["length"] - df_genes["exons_length"]

    print(df_genes[["ID", "num_exons", "length", "exons_length", "introns_length"]])

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate gene statistics")
    parser.add_argument("input", help="Input gff file")
    parser.add_argument("-p", "--prefix", help="Prefix to filter seqid column")
    args = parser.parse_args()
    main(args.input, args.prefix)