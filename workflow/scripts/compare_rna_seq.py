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


rna_df = rna_df.merge(
    genes_df.reset_index()[["index", "ext_gene"]],
    how="left",
    left_on="genes",
    right_on="ext_gene",
).rename(columns={"index": "transcript_index"})
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
