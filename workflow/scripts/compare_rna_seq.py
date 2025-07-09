import pandas as pd
import sys
import altair as alt
import numpy as np

sys.stderr = open(snakemake.log[0], "w")

germ_layer = snakemake.params["germ_layer"]
germ_layer = germ_layer.replace("derm", "")

genes_df = pd.read_csv(snakemake.input["genes_transcripts"], sep="\t")
genes_df = genes_df[genes_df["annotation_type"] == "Promoter"]

rna_df = pd.read_excel(snakemake.input["rna_seq"])  # alternativ z. B. "Sheet1"
rna_df = rna_df[rna_df["test"] == germ_layer].reset_index(drop=True)


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

rna_df = rna_df.groupby("genes", as_index=False).agg(
    {
        "mean_methylation_difference": "mean",
        "absolute_signed_pi_val": "mean",
        "transcript_index": "first",  # oder ggf. "min"/"max", je nach Bedarf
        "ext_gene": "first",  # um eine Spalte für das Gen beizubehalten
        "logFC": "first",  # oder ggf. "median", je nach Bedarf
        "PValue": "first",  # oder ggf. "min"/"max", je nach Bedarf
        "annotation_type": "first",  # um die Annotation beizubehalten
        # weitere Spalten können je nach Bedarf ergänzt werden
    }
)


def plot_scatter(df, title):

    chart = (
        alt.Chart(df)
        .mark_circle()
        .encode(
            x=alt.X("mean_methylation_difference"),
            y=alt.Y("logFC"),
            tooltip=[
                "genes",
            ],
        )
        .properties(
            title=f"{title}",
        )
    )
    return chart


chart_hyper = plot_scatter(
    rna_df[rna_df["mean_methylation_difference"] >= 0], "Hypermethylated"
)
chart_hypo = plot_scatter(
    rna_df[rna_df["mean_methylation_difference"] < 0], "Hypomethylated"
)
chart = (
    alt.hconcat(chart_hyper, chart_hypo)
    .properties(
        title=f"Comparison of RNA-seq logFC to WGS methylation difference for {germ_layer} sites",
    )
    .resolve_scale(y="independent")
)
chart.save(snakemake.output["scatter"])

rna_df[
    [
        "genes",
        "logFC",
        "PValue",
        "transcript_index",
        "annotation_type",
        "mean_methylation_difference",
        "absolute_signed_pi_val",
    ]
].fillna(-2).to_csv(
    snakemake.output["table"],
    sep="\t",
    index=False,
)


#################

# logFC auf eine Nachkommastelle runden
rna_df["logFC_rounded"] = rna_df["logFC"].round(1)
# Gruppen definieren
is_hyper = rna_df["mean_methylation_difference"] >= 0
is_hypo = rna_df["mean_methylation_difference"] < 0
no_dmr = rna_df["mean_methylation_difference"].isna()

# Alle Gene mit Promoter-Annotation

# DataFrame für "nicht in DMRs" (indirekt berechnet)
no_dmr_df = rna_df[no_dmr].copy()
hyper_df = rna_df[is_hyper].copy()
hypo_df = rna_df[is_hypo].copy()


# globales Minimum und Maximum aus gesamten Daten
min_fc = np.floor(rna_df["logFC"].min())
max_fc = np.ceil(rna_df["logFC"].max())

# Bin-Definition: 0.1er Schritte, global festgelegt
bin_settings = alt.Bin(extent=[min_fc, max_fc], step=0.1)


def make_hist(df, title, color):
    return (
        alt.Chart(df)
        .mark_bar(color=color)
        .encode(
            x=alt.X(
                "logFC_rounded:Q",
                bin=bin_settings,  # Globale Bin-Definition
                title="log2FC (rounded to 1 decimal)",
            ),
            y=alt.Y("count()", title="Count"),
            tooltip=[
                alt.Tooltip("logFC_rounded:Q", title="log2FC (rounded to 1 decimal)"),
                alt.Tooltip("count()", title="Count"),
            ],
        )
        .properties(
            title=title,
            width=250,
            height=150,
        )
    )


hist_hyper = make_hist(hyper_df, "Hypermethylated Genes", "#1f77b4")
hist_hypo = make_hist(hypo_df, "Hypomethylated Genes", "#ff7f0e")
hist_none = make_hist(no_dmr_df, "Genes not in DMRs", "#2ca02c")

# Histogramme nebeneinander darstellen
histograms = (
    alt.hconcat(hist_hyper, hist_hypo, hist_none)
    .resolve_scale(x="shared", y="shared")
    .properties(
        title=f"RNA-seq logFC distribution for {germ_layer} sites",
    )
)

# Histogramm speichern
histograms.save(snakemake.output["histogram"])
