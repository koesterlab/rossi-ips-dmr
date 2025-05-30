import pandas as pd
import sys
import altair as alt

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

print(rna_df, file=sys.stderr)


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
chart.save(snakemake.output["plot"])

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
].fillna(-1).to_csv(
    snakemake.output["table"],
    sep="\t",
    index=False,
)
