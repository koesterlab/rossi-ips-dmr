import pandas as pd
import numpy as np

sys.stderr = open(snakemake.log[0], "w", buffering=1)


df = pd.read_csv(snakemake.input[0], sep="\t")


def split_attributes(attributes):
    attr_dict = {}
    for item in attributes.split(";"):
        if "=" in item:
            key, value = item.split("=")
            attr_dict[key] = value
    return attr_dict


attributes_df = df["attributes"].apply(split_attributes).apply(pd.Series)

df = pd.concat([df.drop(columns=["attributes"]), attributes_df], axis=1)


if "gene_name" in df.columns:
    df = df.dropna(subset=["gene_name"])

if "gene_name" in df.columns and "gene_biotype" in df.columns:
    columns = df.columns.tolist()
    gene_name_index = columns.index("gene_name")
    gene_biotype_index = columns.index("gene_biotype")
    columns[gene_name_index], columns[gene_biotype_index] = (
        columns[gene_biotype_index],
        columns[gene_name_index],
    )
    df = df[columns]

df = df[~df["type"].isin(["chromosome", "exon"])]



if "Parent" in df.columns:
    df["Parent"] = df["Parent"].str.replace("gene:", "", regex=True)

df = df[df["q-value"] <= 0.25]

df["absolute_signed_pi_val"] = abs(
    -np.log10(df["p(MWU)"]) * df["mean_methylation_difference"]
)

df = df.sort_values(by="absolute_signed_pi_val", ascending=False)
df = df.dropna(axis=1, how="all")
df.to_csv(snakemake.output[0], sep="\t", index=False)
