import polars as pl
import altair as alt

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
                title="max(qval1, qval2)",
                scale=alt.Scale(range=[30, 1]),
            ),
            tooltip=[
                f"{gene_col}",
                "diffexp",
                # "mean_methylation_difference",
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

    chart.save(output_path)


# DMR and Diffexp information of all genes
comparison_df = pl.read_csv(snakemake.input.comp, separator="\t", null_values="NA")

# Transcription factors with target genes
tf_df = pl.read_csv(snakemake.input.tf_list, separator=",", null_values="NA")

###################### Plot only transcription factors #####################

# The merge ext_gene and source to get all transcription factors of our analysis. We get a new target column, but not all targets are in the comparison df, therefore not every tf has an incfluence. This table gets quite big since every tf can have multiple targets and everz target becomes a new row.
tf_df = comparison_df.join(
    tf_df,
    left_on="ext_gene",
    right_on="source",
    how="inner",
).rename({"ext_gene": "tfs"})


# Keep only rows where the target gene is in the comparison df. This way we lose all targets not included in the diffexp/dmr analysis.
ext_genes = comparison_df.select(pl.col("ext_gene").unique()).to_series().to_list()
tf_df = tf_df.filter(pl.col("target").is_in(ext_genes))
# plot(tf_df, "mean_methylation_difference", "tfs", snakemake.output[0])
tf_df.write_csv(snakemake.output["focus_tfs"], separator="\t")

print(tf_df.filter(pl.col("tfs") == "YY1").filter(pl.col("germ_layer") == "endoderm"))


# We often have multiple DMRs per transcription factor. To sum over all DMRs influencing a target gene, we first need to aggregate the DMRs per tf
tf_df = tf_df.group_by("tfs", "germ_layer", "target", "weight").agg(
    pl.col("mean_methylation_difference").mean(),
    pl.col("qval_dmr").max(),
    pl.col("pval_dmr").max(),
)


# Compute the sum of mean methylation differences of all influencing tfs per target gene and germ layer
tf_effects_per_target = tf_df.group_by("target", "germ_layer").agg(
    (pl.col("weight") * pl.col("mean_methylation_difference"))
    # (pl.col("mean_methylation_difference"))
    .sum().alias("tf_sum_mean_methylation_difference"),
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


###################### Prepare for datavzrd #####################
comparison_with_tf = (
    comparison_with_tf.with_columns(
        (pl.col("mean_methylation_difference_tf_adjusted") * pl.col("diffexp")).alias(
            "methylation_diff_x_diffexp"
        )
    )
    .with_columns(
        pl.when(
            (
                pl.col("mean_methylation_difference_tf_adjusted").sign()
                == pl.col("diffexp").sign()
            )
        )
        .then(pl.col("methylation_diff_x_diffexp") * 0.5)
        .when(pl.col("mean_methylation_difference_tf_adjusted").sign() >= 0)
        .then(pl.col("methylation_diff_x_diffexp") * -1)
        .otherwise(pl.col("methylation_diff_x_diffexp"))
        .alias("ranked_meth_diffexp")
    )
    .sort("methylation_diff_x_diffexp", descending=False)
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
            "methylation_diff_x_diffexp",
            "ranked_meth_diffexp",
        ]
    )
    .rename({"mean_methylation_difference": "mean_methylation_difference_original"})
    .with_row_count("row_id")
)

plot(
    comparison_with_tf,
    "mean_methylation_difference_tf_adjusted",
    "ext_gene",
    snakemake.output["plot"],
)
# Are there any inf or NA values in the qval_combined column?
comparison_with_tf = comparison_with_tf.with_columns(
    pl.col("qval_combined").fill_null(1.0).fill_nan(1.0)
)
comparison_with_tf.write_csv(snakemake.output["comp_tf_adj"], separator="\t")
