import polars as pl
import altair as alt

pl.Config.set_tbl_rows(10)
pl.Config.set_tbl_cols(30)


filename_to_name = {
    "distal_intergenic": "Distal Intergenic",
    "promoter": "Promoter",
    "intron": "Intron",
    "exon": "Exon",
    "3_utr": "3' UTR",
    "5_utr": "5' UTR",
    "downstream": "Downstream",
}
annotation_type = filename_to_name.get(snakemake.params.get("annotation_type", None))


def read_dmrs(path: str, layer: str) -> pl.DataFrame:
    return pl.read_csv(
        path,
        separator="\t",
        schema_overrides={"chr": pl.Utf8},
        null_values="NA",
    ).with_columns(pl.lit(layer).alias("germ_layer"))


diffexp = snakemake.input.diffexp
ectoderm = snakemake.input.ectoderm
endoderm = snakemake.input.endoderm
mesoderm = snakemake.input.mesoderm


diffexp_df = pl.read_csv(diffexp, separator="\t", null_values="NA").rename(
    {
        "qval": "qval_diffexp",
    }
)

ectoderm_df = read_dmrs(ectoderm, "ectoderm")
endoderm_df = read_dmrs(endoderm, "endoderm")
mesoderm_df = read_dmrs(mesoderm, "mesoderm")
dmrs_df = (
    pl.concat(
        [
            ectoderm_df,
            endoderm_df,
            mesoderm_df,
        ]
    )
    .unique()
    .rename(
        {
            "qval": "qval_dmr",
        }
    )
)

common_df = diffexp_df.join(
    dmrs_df,
    on="ext_gene",
    how="inner",
)
common_df = common_df.filter(pl.col("annotation_type") == annotation_type).with_columns(
    (pl.max_horizontal(pl.col("qval_diffexp"), pl.col("qval_dmr"))).alias(
        "qval_combined"
    )
)
common_df = common_df.with_columns(
    pl.when(pl.col("germ_layer") == "ectoderm")
    .then(pl.col("b_conditionectoderm"))
    .when(pl.col("germ_layer") == "endoderm")
    .then(pl.col("b_conditionendoderm"))
    .when(pl.col("germ_layer") == "mesoderm")
    .then(pl.col("b_conditionmesoderm"))
    .otherwise(None)
    .alias("diffexp")
)

layer_select = alt.selection_point(
    fields=["germ_layer"], bind="legend", name="Germ Layer"
)
color = alt.condition(
    layer_select,
    alt.Color(
        "germ_layer:N",
        title="Germ Layer",
        scale=alt.Scale(scheme="category10"),
    ),
    alt.value("lightgray"),
)


chart = (
    alt.Chart(common_df.to_pandas())
    .mark_point(filled=True)
    .encode(
        x="diffexp:Q",
        y="mean_methylation_difference:Q",
        size=alt.Size(
            "qval_combined:Q",
            title="max(qval1, qval2)",
            scale=alt.Scale(range=[30, 1]),
        ),
        tooltip=[
            "germ_layer",
            "ens_gene",
            "diffexp",
            "mean_methylation_difference",
            "qval_combined",
            "qval_diffexp",
            "qval_dmr",
        ],
        color=color,
        opacity=alt.Opacity(
            "qval_combined:Q",
            scale=alt.Scale(range=[1, 0.1]),
            title="max(qval1, qval2)",
        ),
    )
    .add_params(layer_select)
)

chart.save(snakemake.output[0])
