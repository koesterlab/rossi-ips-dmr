
import polars as pl

import altair as alt
import polars as pl

alt.data_transformers.disable_max_rows()

pl.Config.set_tbl_rows(10)
pl.Config.set_tbl_cols(100)

ANNOTATION_TYPE_NAMES = {
    "distal_intergenic": "Distal Intergenic",
    "promoter": "Promoter",
    "intron": "Intron",
    "exon": "Exon",
    "3_utr": "3' UTR",
    "5_utr": "5' UTR",
    "downstream": "Downstream",
}

annotation_type = ANNOTATION_TYPE_NAMES.get(snakemake.params.get("annotation"), None)
meth_df = pl.read_csv(
    snakemake.input.meth,
    separator="\t",
    null_values=["NA"],
    infer_schema_length=None,  # scans entire file to infer types
).select("chromosome", "position", "annotation", "transcriptId", "psc_methylation", "endoderm_methylation", "ectoderm_methylation", "mesoderm_methylation").with_columns(
    pl.col("annotation").str.split(" ").list.get(0)
)
expr_df = pl.read_csv(
    snakemake.input.expr,
    separator="\t",
    null_values=["NA"],
    infer_schema_length=None,
).with_columns(
    pl.col("transcript").str.split(".").list.get(0).alias("transcriptId")
)
print(annotation_type)
if annotation_type != None:
    meth_df = meth_df.filter(pl.col("annotation") == annotation_type)


expr_df = expr_df.with_columns(
    pl.mean_horizontal(pl.col("^ectoderm.*$")).alias("ectoderm_expression"),
    pl.mean_horizontal(pl.col("^endoderm.*$")).alias("endoderm_expression"),
    pl.mean_horizontal(pl.col("^mesoderm.*$")).alias("mesoderm_expression"),
    pl.mean_horizontal(pl.col("^psc.*$")).alias("psc_expression"),
).select(
    "transcriptId",
    "gene",
    "ectoderm_expression",
    "endoderm_expression",
    "mesoderm_expression",
    "psc_expression",
)
meth_df = (
    meth_df
    .group_by("transcriptId")
    .agg(
        pl.col("psc_methylation").mean(),
        pl.col("endoderm_methylation").mean(),
        pl.col("ectoderm_methylation").mean(),
        pl.col("mesoderm_methylation").mean(),
    )
)


df = meth_df.join(expr_df, on="transcriptId")



LAYERS = [
    ("psc_expression", "psc_methylation", "PSC"),
    ("endoderm_expression", "endoderm_methylation", "Endoderm"),
    ("ectoderm_expression", "ectoderm_methylation", "Ektoderm"),
    ("mesoderm_expression", "mesoderm_methylation", "Mesoderm"),
]

long_df = pl.concat(
    [
        df.select(
            pl.col("transcriptId"),
            pl.col("gene"),
            pl.col(expr_col).alias("expression"),
            pl.col(meth_col).alias("methylation"),
            pl.lit(label).alias("layer"),
        ).drop_nulls(["expression", "methylation"])
        for expr_col, meth_col, label in LAYERS
    ]
).with_columns(
    (pl.col("expression") + 1).log(base=2).alias("log2_expression")
)


chart = (
    alt.Chart(long_df.to_pandas())
    .mark_circle(size=25, opacity=0.4)
    .encode(
        x=alt.X("log2_expression:Q", title="log₂(Expression + 1)"),
        y=alt.Y("methylation:Q", title="Methylierung (%)"),
        tooltip=["gene:N", "transcriptId:N", "expression:Q", "methylation:Q"],
        facet=alt.Facet(
            "layer:N",
            title=None,
            header=alt.Header(labelFontSize=13, labelFontWeight="bold"),
        ),
    )
    .properties(width=250, height=250)
    .resolve_scale(x="independent", y="independent")
    .properties(
        title=alt.TitleParams(
            "Expression vs. Methylierung pro Keimblatt",
            fontSize=16,
            fontWeight="bold",
        )
    )
)

chart.save(snakemake.output[0])
