import sys

import altair as alt
import polars as pl

sys.stderr = open(snakemake.log[0], "w", buffering=1)

pl.Config.set_tbl_rows(10)
pl.Config.set_tbl_cols(300)


def plot(df, effect_col, gene_col, output_path):
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
        alt.Chart(df.to_pandas())
        .mark_point(filled=True)
        .encode(
            x="diffexp:Q",
            y=f"{effect_col}:Q",
            size=alt.Size(
                "qval_combined:Q",
                title="min(qDMR, qDGE)",
                scale=alt.Scale(range=[30, 1]),
            ),
            tooltip=[
                f"{gene_col}",
                "diffexp",
                "qval_combined",
                "qval_diffexp",
                "qval_dmr",
            ],
            color=color,
            opacity=alt.Opacity(
                "qval_combined:Q",
                scale=alt.Scale(range=[1, 0.1]),
                title="min(qDMR, qDGE)",
            ),
        )
        .add_params(layer_select)
    )

    chart.save(output_path)


# DMR and Diffexp information of all genes
comparison_df = pl.read_csv(snakemake.input.comp, separator="\t", null_values="NA")

# Transcription factors with target genes
tf_df = pl.read_csv(snakemake.input.tf_list, separator=",", null_values="NA")

###################### Plot only transcription factors #####################

# Join on the TF name to get all TFs of our analysis with their targets. Every
# TF-target pair becomes a separate row, so this table gets quite big.
tf_df = comparison_df.join(
    tf_df,
    left_on="ext_gene",
    right_on="source",
    how="inner",
).rename({"ext_gene": "tfs"})

# Keep only targets that are part of the diffexp/DMR analysis.
tf_df = tf_df.filter(pl.col("target").is_in(comparison_df["ext_gene"].unique()))


# A TF often has multiple DMRs; aggregate them per TF before summing over all
# TFs influencing a target gene.
tf_df = tf_df.group_by("tfs", "germ_layer", "target", "weight").agg(
    pl.col("mean_methylation_difference").mean(),
    pl.col("qval_dmr").max(),
    pl.col("pval_dmr").max(),
)
tf_df.write_csv(snakemake.output["focus_tfs"], separator="\t")


# Compute the sum of mean methylation differences of all influencing tfs per target gene and germ layer
tf_effects_per_target = tf_df.group_by("target", "germ_layer").agg(
    (pl.col("weight") * pl.col("mean_methylation_difference"))
    .sum()
    .alias("tf_sum_mean_methylation_difference"),
    pl.col("tfs").unique().sort().str.join(","),
)


# Merge comparison df with TF target info
comparison_with_tf = comparison_df.join(
    tf_effects_per_target,
    left_on=["ext_gene", "germ_layer"],
    right_on=["target", "germ_layer"],
    how="left",
)

# Choose the sum if the target is influenced by tfs, otherwise keep the original mean methylation difference
comparison_with_tf = comparison_with_tf.with_columns(
    pl.when(pl.col("tf_sum_mean_methylation_difference").is_not_null())
    .then(pl.col("tf_sum_mean_methylation_difference"))
    .otherwise(pl.col("mean_methylation_difference"))
    .alias("mean_methylation_difference_tf_adjusted")
).with_columns(
    pl.col("tf_sum_mean_methylation_difference").is_not_null().alias("is_tf_target")
)

# Symmetric axis ranges around 0; methylation differences are scaled to the
# diffexp range so that both axes contribute equally to the distances below.
diffexp_max, meth_diff_max = comparison_with_tf.select(
    pl.col("diffexp").abs().max(),
    pl.col("mean_methylation_difference_tf_adjusted").abs().max(),
).row(0)
diffexp_min = -diffexp_max
meth_diff_max_scaled = meth_diff_max * diffexp_max
meth_diff_min_scaled = -meth_diff_max_scaled

max_dist = ((diffexp_max) ** 2 + meth_diff_max_scaled**2) ** 0.5


###################### Prepare for datavzrd #####################
comparison_with_tf = (
    comparison_with_tf.with_columns(
        (pl.col("mean_methylation_difference_tf_adjusted") * diffexp_max).alias(
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
            "mean_methylation_difference_tf_adjusted",
            "tfs",
            "is_tf_target",
            "ranked_meth_diffexp",
        ]
    )
    .with_row_index("row_id")
)

plot(
    comparison_with_tf,
    "mean_methylation_difference_tf_adjusted",
    "ext_gene",
    snakemake.output["plot"],
)
# Missing q-values (gene only present in one of both analyses) count as not significant
comparison_with_tf = comparison_with_tf.with_columns(
    pl.col("qval_combined").fill_null(1.0).fill_nan(1.0)
)
comparison_with_tf.write_csv(snakemake.output["comp_tf_adj"], separator="\t")
