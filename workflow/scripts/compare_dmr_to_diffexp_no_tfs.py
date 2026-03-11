import sys

import altair as alt
import polars as pl

sys.stderr = open(snakemake.log[0], "w", buffering=1)

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

# The three non-base layers are passed in as a list from the rule, e.g.
# ["ectoderm", "endoderm", "mesoderm"] when base == "psc".
non_base_layers: list[str] = snakemake.params["non_base_layers"]

# The input keys layer1/layer2/layer3 map to the three non-base layers
# in the same order as non_base_layers.
layer_inputs = [snakemake.input.layer1, snakemake.input.layer2, snakemake.input.layer3]


def read_dmrs(path: str, layer: str) -> pl.DataFrame:
    return pl.read_csv(
        path,
        separator="\t",
        schema_overrides={"chr": pl.Utf8},
        null_values="NA",
    ).with_columns(pl.lit(layer).alias("germ_layer"))


diffexp_df = pl.read_csv(
    snakemake.input.diffexp, separator="\t", null_values="NA"
).rename(
    {
        "qval": "qval_diffexp",
        "pval": "pval_diffexp",
    }
)

dmrs_df = (
    pl.concat(
        [read_dmrs(path, layer) for path, layer in zip(layer_inputs, non_base_layers)]
    )
    .unique()
    .rename(
        {
            "qval": "qval_dmr",
            "pval": "pval_dmr",
        }
    )
)

common_df = diffexp_df.join(dmrs_df, on="ext_gene", how="inner").drop_nulls(
    [
        "mean_methylation_difference",
        "qval_diffexp",
        "qval_dmr",
        "pval_diffexp",
        "pval_dmr",
    ]
)

common_df = (
    common_df.filter(pl.col("annotation_type") == annotation_type)
    .with_columns(
        (pl.max_horizontal(pl.col("qval_diffexp"), pl.col("qval_dmr"))).alias(
            "qval_combined"
        )
    )
    .with_columns(
        (pl.max_horizontal(pl.col("pval_diffexp"), pl.col("pval_dmr"))).alias(
            "pval_combined"
        )
    )
)

# Build the diffexp / diffexp_se columns dynamically based on which layers
# are actually present.  Sleuth names the beta-value column
# b_condition<layer> (relative to the chosen base_level).
beta_expr = pl.lit(None, dtype=pl.Float64)
beta_se_expr = pl.lit(None, dtype=pl.Float64)

for layer in non_base_layers:
    beta_col = f"b_condition{layer}"
    beta_se_col = f"b_condition{layer}_se"
    beta_expr = (
        pl.when(pl.col("germ_layer") == layer)
        .then(pl.col(beta_col))
        .otherwise(beta_expr)
    )
    beta_se_expr = (
        pl.when(pl.col("germ_layer") == layer)
        .then(pl.col(beta_se_col))
        .otherwise(beta_se_expr)
    )

common_df = common_df.with_columns(
    beta_expr.alias("diffexp"),
    beta_se_expr.alias("diffexp_se"),
)
common_df = common_df.filter(pl.col("diffexp").is_not_null())

print(common_df.head(6))
common_df = common_df.group_by(
    [c for c in common_df.columns if c not in ["ens_gene", "target_id", "mane"]]
).agg(pl.col("ens_gene").first().alias("ens_gene"))
print(common_df.head(6))

x_domain = [
    common_df["diffexp"].min(),
    common_df["diffexp"].max(),
]

y_domain = [
    common_df["mean_methylation_difference"].min(),
    common_df["mean_methylation_difference"].max(),
]


layer_select = alt.selection_point(
    fields=["germ_layer"], bind="legend", name="Germ Layer"
)

qval_slider = alt.param(
    name="qval_max",
    value=0.5,
    bind=alt.binding_range(min=0, max=1, step=0.01, name="max qval: "),
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
    .transform_filter(alt.datum.qval_combined <= qval_slider)
    .mark_point(filled=True)
    .encode(
        x=alt.X("diffexp:Q", scale=alt.Scale(domain=x_domain)),
        y=alt.Y("mean_methylation_difference:Q", scale=alt.Scale(domain=y_domain)),
        size=alt.Size(
            "qval_combined:Q",
            title="max(qval1, qval2)",
            scale=alt.Scale(range=[30, 1]),
        ),
        tooltip=[
            "ext_gene",
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
    .add_params(layer_select, qval_slider)
)

chart.save(snakemake.output.html)

diffexp_min, diffexp_max, meth_diff_min, meth_diff_max = common_df.select(
    pl.col("diffexp").min().alias("x_min"),
    pl.col("diffexp").max().alias("x_max"),
    pl.col("mean_methylation_difference").min().alias("y_min"),
    pl.col("mean_methylation_difference").max().alias("y_max"),
).row(0)

(diffexp_min, diffexp_max) = (
    -max(abs(diffexp_min), abs(diffexp_max)),
    max(abs(diffexp_min), abs(diffexp_max)),
)
(meth_diff_min, meth_diff_max) = (
    -max(abs(meth_diff_min), abs(meth_diff_max)),
    max(abs(meth_diff_min), abs(meth_diff_max)),
)
meth_diff_max_scaled = meth_diff_max * diffexp_max
meth_diff_min_scaled = meth_diff_min * diffexp_max

max_dist = ((diffexp_max) ** 2 + meth_diff_max_scaled**2) ** 0.5

common_df = (
    common_df.with_columns(
        (pl.col("mean_methylation_difference") * diffexp_max).alias(
            "mean_methylation_difference_scaled"
        )
    )
    .with_columns(
        (
            (
                (pl.col("diffexp") - diffexp_min) ** 2
                + (pl.col("mean_methylation_difference_scaled") - meth_diff_max_scaled)
                ** 2
            ).sqrt()
        ).alias("dist_top_left"),
        (
            (
                (pl.col("diffexp") - diffexp_max) ** 2
                + (pl.col("mean_methylation_difference_scaled") - meth_diff_min_scaled)
                ** 2
            ).sqrt()
        ).alias("dist_bottom_right"),
    )
    .with_columns(
        pl.when(pl.col("dist_top_left") < pl.col("dist_bottom_right"))
        .then(max_dist - pl.col("dist_top_left"))
        .otherwise(-max_dist + pl.col("dist_bottom_right"))
        .alias("ranked_meth_diffexp")
    )
    .sort(pl.col("ranked_meth_diffexp").abs(), descending=True)
    .select(
        [
            "ext_gene",
            "ens_gene",
            "germ_layer",
            "qval_diffexp",
            "qval_dmr",
            "qval_combined",
            "pval_dmr",
            "pval_diffexp",
            "diffexp",
            "diffexp_se",
            "mean_methylation_difference",
            "ranked_meth_diffexp",
        ]
    )
    .with_row_index("row_id")
)


common_df.write_csv(snakemake.output.tsv, separator="\t")
