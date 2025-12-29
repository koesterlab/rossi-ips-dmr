import polars as pl
import altair as alt

pl.Config.set_tbl_rows(10)
pl.Config.set_tbl_cols(30)


def read_dmrs(path: str, layer: str) -> pl.DataFrame:
    return (
        pl.read_csv(
            path,
            separator="\t",
            schema_overrides={"chr": pl.Utf8},
        )
        .filter(pl.col("annotation_type") == "Promoter")
        .with_columns(pl.lit(layer).alias("germ_layer"))
    )


diffexp = snakemake.input.diffexp
ectoderm = snakemake.input.ectoderm
endoderm = snakemake.input.endoderm
mesoderm = snakemake.input.mesoderm


diffexp_df = pl.read_csv(diffexp, separator="\t").rename(
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
    left_on="ext_gene",
    right_on="ext_gene",
    how="inner",
)

common_df = common_df.with_columns(
    (1 - pl.col("qval_diffexp") * pl.col("qval_dmr")).alias("qval_combined")
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
print(common_df.head())

layer_select = alt.selection_point(fields=["germ_layer"], bind="legend")


chart = (
    alt.Chart(common_df.to_pandas())
    .mark_point()
    .encode(
        x="diffexp:Q",
        y="mean_methylation_difference:Q",
        tooltip=[
            "germ_layer",
            "ens_gene",
            "diffexp",
            "mean_methylation_difference",
            "qval_combined",
            "qval_diffexp",
            "qval_dmr",
        ],
        color=alt.Color(
            "germ_layer:N",
            title="Germ Layer",
            scale=alt.Scale(scheme="category10"),
        ),
        opacity=alt.condition(
            layer_select,
            alt.Opacity(
                "qval_combined:Q",
                scale=alt.Scale(domain=[0.8, 1]),
                title="1 − combined q-value",
            ),
            alt.OpacityValue(0.05),
        ),
    )
    .add_params(layer_select)
)


chart.save(snakemake.output[0])
