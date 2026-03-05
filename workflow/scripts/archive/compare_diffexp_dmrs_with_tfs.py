import polars as pl
import altair as alt

# Polars Anzeigeoptionen
pl.Config.set_tbl_rows(500)
pl.Config.set_tbl_cols(300)

# Mapping für Annotation Types
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


# Read DMRs
def read_dmrs(path: str, layer: str) -> pl.DataFrame:
    return (
        pl.read_csv(
            path,
            separator="\t",
            schema_overrides={"chr": pl.Utf8},
            null_values="NA",
        )
        .with_columns(pl.lit(layer).alias("germ_layer"))
        .filter(pl.col("annotation_type") == annotation_type)
        .group_by(["ext_gene", "germ_layer"])
        .agg(
            [
                pl.mean("mean_methylation_difference"),
                pl.mean("qval"),
                pl.mean("absolute_signed_pi_val"),
            ]
        )
        .select(
            [
                "ext_gene",
                "mean_methylation_difference",
                "qval",
                "germ_layer",
                "absolute_signed_pi_val",
            ]
        )
    )


# Differential-Expression Daten
diffexp_df = (
    pl.read_csv(snakemake.input.diffexp, separator="\t", null_values="NA")
    .rename({"qval": "qval_diffexp"})
    .select(
        [
            "ext_gene",
            "ens_gene",
            "b_conditionectoderm",
            "b_conditionectoderm_se",
            "b_conditionendoderm",
            "b_conditionendoderm_se",
            "b_conditionmesoderm",
            "b_conditionmesoderm_se",
            "qval_diffexp",
        ]
    )
)


# --- TF DMR Aggregation ---
def tf_dmr_sum(tf_path: str, dmrs_df: pl.DataFrame, layer: str) -> pl.DataFrame:
    tf_df = (
        pl.read_csv(tf_path, separator=",")
        .rename(
            {
                "TFs": "TF_list",
                "target": "ext_gene",
            }
        )
        .with_columns(pl.col("TF_list").str.split(", ").alias("TF"))
        .explode("TF")  # put each TF in its own row
        .join(
            dmrs_df,  # join with dmr dataset to get methylation differences per ext_gene
            on="ext_gene",
            how="right",
        )
    )
    # Join another time to get methylation differences per TF. The TF corresponds to the original ext_gene. We do not find a value for each transcription factor, since the TFs are from the diffexp and we merge with DMRs here.
    tf_dmr = tf_df.join(
        dmrs_df.select(
            pl.col("ext_gene").alias("TF"),
            pl.col("mean_methylation_difference").alias(
                "mean_methylation_difference_TF"
            ),
        ),
        on="TF",
        how="left",
    )

    # Now aggregate per target gene. The mean_methylation_difference_TF is the summed value of all DMRs associated with TFs targeting that gene.
    return tf_dmr.group_by("ext_gene").agg(
        [
            pl.sum("mean_methylation_difference_TF").alias("dmr_sum_TFs"),
            pl.first("qval").alias("qval"),
            pl.first("absolute_signed_pi_val").alias("mean_absolute_signed_pi_val"),
            pl.first("germ_layer").alias("germ_layer"),
            pl.first("mean_methylation_difference"),
        ]
    )


# DMR files per germ layer
ectoderm_df = read_dmrs(snakemake.input.ectoderm, "ectoderm")
endoderm_df = read_dmrs(snakemake.input.endoderm, "endoderm")
mesoderm_df = read_dmrs(snakemake.input.mesoderm, "mesoderm")

# taget genes with their TFs
ectoderm_tfs = snakemake.input.ectoderm_tfs
endoderm_tfs = snakemake.input.endoderm_tfs
mesoderm_tfs = snakemake.input.mesoderm_tfs

ectoderm_dmr_sum_from_tfs = tf_dmr_sum(ectoderm_tfs, ectoderm_df, "ectoderm")
endoderm_dmr_sum_from_tfs = tf_dmr_sum(endoderm_tfs, endoderm_df, "endoderm")
mesoderm_dmr_sum_from_tfs = tf_dmr_sum(mesoderm_tfs, mesoderm_df, "mesoderm")


# Concat the DMRs from all germ layers
dmrs_df = (
    pl.concat(
        [
            ectoderm_dmr_sum_from_tfs,
            endoderm_dmr_sum_from_tfs,
            mesoderm_dmr_sum_from_tfs,
        ]
    )
    .unique()
    .rename({"qval": "qval_dmr"})
    .select(
        [
            "ext_gene",
            "mean_methylation_difference",
            "dmr_sum_TFs",
            "qval_dmr",
            "germ_layer",
        ]
    )
    .group_by("ext_gene", "germ_layer")
    .agg(
        [
            pl.mean("mean_methylation_difference"),
            pl.mean("qval_dmr"),
            pl.sum("dmr_sum_TFs"),
        ]
    )
)


common_df = diffexp_df.join(dmrs_df, on="ext_gene", how="inner").filter(
    (
        pl.col("mean_methylation_difference").is_not_null()
        | pl.col("dmr_sum_TFs").is_not_null()
    )
    & pl.col("qval_diffexp").is_not_null()
    & pl.col("qval_dmr").is_not_null()
)

# Combine q-values
common_df = common_df.with_columns(
    (pl.max_horizontal(pl.col("qval_diffexp"), pl.col("qval_dmr"))).alias(
        "qval_combined"
    )
)

# Choose diffexp columns based on germ layer
common_df = (
    common_df.with_columns(
        pl.when(pl.col("germ_layer") == "ectoderm")
        .then(pl.col("b_conditionectoderm"))
        .when(pl.col("germ_layer") == "endoderm")
        .then(pl.col("b_conditionendoderm"))
        .when(pl.col("germ_layer") == "mesoderm")
        .then(pl.col("b_conditionmesoderm"))
        .otherwise(None)
        .alias("diffexp")
    )
    .with_columns(
        pl.when(pl.col("germ_layer") == "ectoderm")
        .then(pl.col("b_conditionectoderm_se"))
        .when(pl.col("germ_layer") == "endoderm")
        .then(pl.col("b_conditionendoderm_se"))
        .when(pl.col("germ_layer") == "mesoderm")
        .then(pl.col("b_conditionmesoderm_se"))
        .otherwise(None)
        .alias("diffexp_se")
    )
    .filter(pl.col("diffexp").is_not_null())
)

# If we have a DMR sum from TFs, use that as methylation difference, otherwise use the mean methylation difference of the gene itself
common_df = common_df.with_columns(
    pl.when(pl.col("dmr_sum_TFs").is_not_null() & (pl.col("dmr_sum_TFs") != 0))
    .then(pl.col("dmr_sum_TFs"))
    .otherwise(pl.col("mean_methylation_difference"))
    .alias("meth_difference_per_gene")
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
        y="meth_difference_per_gene:Q",
        size=alt.Size(
            "qval_combined:Q",
            title="max(qval1, qval2)",
            scale=alt.Scale(range=[30, 1]),
        ),
        tooltip=[
            "ext_gene",
            "diffexp",
            "meth_difference_per_gene",
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


# --- Optional: export table ---
common_df = (
    # common_df.filter(pl.col("qval_combined") <= 0.5)
    common_df.with_columns(
        (pl.col("meth_difference_per_gene") * pl.col("diffexp")).alias(
            "dmr_sum_x_diffexp"
        )
    )
    .sort("dmr_sum_x_diffexp", descending=False)
    .select(
        [
            "ext_gene",
            "ens_gene",
            "germ_layer",
            "qval_diffexp",
            "qval_dmr",
            "diffexp",
            "diffexp_se",
            "mean_methylation_difference",
            # "methylation_diff_x_diffexp",
            "dmr_sum_x_diffexp",
        ]
    )
)


common_df.write_csv(snakemake.output[1], separator="\t")
