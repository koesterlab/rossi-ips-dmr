import pandas as pd
import sys

sys.stderr = open(snakemake.log[0], "w")

# Pfad zur Excel-Datei
germ_layer = snakemake.params["germ_layer"]

genes_df = pd.read_csv(snakemake.input["genes_transcripts"], sep="\t")

# Alle verfügbaren Tabellenblätter anzeigen
# Ein bestimmtes Tabellenblatt einlesen (z. B. erstes Sheet)
rna_df = pd.read_excel(snakemake.input["rna_seq"])  # alternativ z. B. "Sheet1"
print(rna_df.head(), file=sys.stderr)
germ_layer = germ_layer.replace("derm", "")
print(germ_layer, file=sys.stderr)

# When germ_layer is "endoderm" make endo out of it

rna_df = rna_df[rna_df["test"] == germ_layer].reset_index(drop=True)
# Zeige die ersten 5 Zeilen zur Kontrolle
print(genes_df.head(), file=sys.stderr)
print(rna_df.head(), file=sys.stderr)


gene_to_index = (
    genes_df.reset_index()
    .drop_duplicates(subset="ext_gene")
    .set_index("ext_gene")["index"]
)
rna_df["transcript_index"] = rna_df["genes"].map(gene_to_index).astype("Int64")

rna_df = pd.merge(
    rna_df,
    genes_df,
    left_on="genes",
    right_on="ext_gene",
    how="left",
)
# Ausgabe zur Kontrolle
print(rna_df[["genes", "transcript_index"]].head())


rna_df[
    [
        "genes",
        "PValue",
        "transcript_index",
        "annotation_type",
        "mean_methylation_difference",
        "absolute_signed_pi_val",
    ]
].to_csv(
    snakemake.output[0],
    sep="\t",
    index=False,
)
