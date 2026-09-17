import sys

import altair as alt
import polars as pl

sys.stderr = open(snakemake.log[0], "w", buffering=1)
alt.data_transformers.disable_max_rows()
alt.data_transformers.enable("vegafusion")

LAYERS = ["psc", "endoderm", "ectoderm", "mesoderm"]
ANNOTATION_TYPE_NAMES = {
    "distal_intergenic": "Distal Intergenic",
    "promoter": "Promoter",
    "intron": "Intron",
    "exon": "Exon",
    "3_utr": "3' UTR",
    "5_utr": "5' UTR",
    "downstream": "Downstream",
}
# Any other value (e.g. "all") disables the annotation filter
annotation_type = ANNOTATION_TYPE_NAMES.get(snakemake.wildcards.annotation)

meth_df = pl.read_csv(
    snakemake.input.meth,
    separator="\t",
    null_values=["NA"],
    columns=["annotation", "transcriptId", *[f"{layer}_methylation" for layer in LAYERS]],
    infer_schema_length=None,
).with_columns(
    # e.g. "Promoter (<=1kb)" -> "Promoter", "Distal Intergenic" stays as is
    pl.col("annotation").str.replace(r"\s*\(.*\)$", "")
)
if annotation_type is not None:
    meth_df = meth_df.filter(pl.col("annotation") == annotation_type)

# Mean methylation of all CpGs per transcript region.
meth_df = meth_df.group_by("transcriptId", "annotation").agg(
    pl.col("psc_methylation").mean(),
    pl.col("endoderm_methylation").mean(),
    pl.col("ectoderm_methylation").mean(),
    pl.col("mesoderm_methylation").mean(),
)

# Mean expression over all replicates (columns are named "<layer>_<replicate>")
expr_df = pl.read_csv(
    snakemake.input.expr, separator="\t", null_values=["NA"], infer_schema_length=None
).select(
    pl.col("transcript").str.split(".").list.get(0).alias("transcriptId"),
    *[pl.mean_horizontal(pl.col(f"^{layer}_.*$")).alias(f"{layer}_expression") for layer in LAYERS],
)


df = meth_df.join(expr_df, on="transcriptId")
long_df = (
    pl.concat(
        df.select(
            pl.lit(layer).alias("layer"),
            pl.col(f"{layer}_expression").alias("expression"),
            pl.col(f"{layer}_methylation").alias("methylation"),
        )
        for layer in LAYERS
    )
    .drop_nulls()
    .with_columns(log2_expression=(pl.col("expression") + 1).log(base=2))
)

chart = (
    alt.Chart(long_df.to_pandas())
    .mark_rect()
    .encode(
        x=alt.X("log2_expression:Q", bin=alt.Bin(maxbins=100), title="log₂(TPM + 1)"),
        y=alt.Y("methylation:Q", bin=alt.Bin(maxbins=100), title="Methylation (%)"),
        color=alt.Color(
            "count():Q",
            scale=alt.Scale(type="log"),
            title="Count",
            legend=alt.Legend(format=",d"),
        ),
        facet=alt.Facet(
            "layer:N",
            title=None,
            header=alt.Header(labelFontSize=13, labelFontWeight="bold"),
        ),
    )
    .properties(width=250, height=250)
)
chart.save(snakemake.output[0])
