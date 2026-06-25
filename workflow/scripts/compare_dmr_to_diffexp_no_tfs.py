import sys

import altair as alt
import polars as pl

sys.stderr = open(snakemake.log[0], "w", buffering=1)

pl.Config.set_tbl_rows(100)
pl.Config.set_tbl_cols(100)

LAYER_COLORS = {
    "endoderm": "#ff7f0e",
    "mesoderm": "#2ca02c",
    "ectoderm": "#d62728",
}

ANNOTATION_TYPE_NAMES = {
    "distal_intergenic": "Distal Intergenic",
    "promoter": "Promoter",
    "intron": "Intron",
    "exon": "Exon",
    "3_utr": "3' UTR",
    "5_utr": "5' UTR",
    "downstream": "Downstream",
}

annotation_type = ANNOTATION_TYPE_NAMES.get(snakemake.params.get("annotation_type"))
non_base_layers: list[str] = snakemake.params["non_base_layers"]
diffexp_inputs = snakemake.input.diffexp
# Direction of diffexp depends on the order of conditions in the model name.
# If the requested direction doesn't match how the model is named in the
# config, the table is still returned as defined in the config, but its
# diffexp values need to be multiplied by -1 to match the requested direction.
diffexp_base_signs = snakemake.params["diffexp_base_signs"]
layer_inputs = [snakemake.input.layer1, snakemake.input.layer2, snakemake.input.layer3]

DMR_COLUMNS = [
    "ext_gene",
    "transcriptId",
    "ens_gene",
    "annotation_type",
    "mean_methylation_difference",
    "absolute_signed_pi_val",
    "qval_dmr",
    "pval_dmr",
    "chr",
    "start_dmr",
    "end_dmr",
    "germ_layer",
]

DIFFEXP_COLUMNS = [
    "ext_gene",
    "ens_gene",
    "target_id",
    "qval_diffexp",
    "pval_diffexp",
    "diffexp",
    "diffexp_se",
    "signed_pi_col",
    "germ_layer",
]


def read_dmrs(path: str, layer: str) -> pl.DataFrame:
    """Load one germ layer's DMR table and tag it with its layer name."""
    return (
        pl.read_csv(
            path,
            separator="\t",
            schema_overrides={"chr": pl.Utf8},
            null_values="NA",
        )
        .rename({"qval": "qval_dmr", "pval": "pval_dmr"})
        .with_columns(germ_layer=pl.lit(layer))
        # Not filtered by annotation_type here: find_val_genes() needs the
        # unfiltered table to check candidate genes against all annotations.
        .select(DMR_COLUMNS)
    )


def read_diffexp(path: str, layer: str, sign: int) -> pl.DataFrame:
    """
    Load one germ layer's differential expression table, tag it with its
    layer name, and flip the sign of its effect-size columns if `sign == -1`
    (see diffexp_base_signs comment above).
    """
    df = pl.read_csv(path, separator="\t", null_values="NA")

    # There is only one column mathing for each column under interest.
    effect_col = next(
        c for c in df.columns if c.startswith("b_") and c.endswith(("-", "+"))
    )
    signed_pi_col = next(
        c for c in df.columns if c.startswith("signed_pi") and c.endswith(("-", "+"))
    )
    se_col = next(c for c in df.columns if c.endswith("_se"))

    return (
        df.rename({"qval": "qval_diffexp", "pval": "pval_diffexp"})
        .with_columns(
            diffexp=pl.col(effect_col) * sign,
            signed_pi_col=pl.col(signed_pi_col) * sign,
            diffexp_se=pl.col(se_col),
            germ_layer=pl.lit(layer),
        )
        .select(DIFFEXP_COLUMNS)
    )


def write_val_genes_table(
    combined_df: pl.DataFrame,
) -> None:
    """Write the per-validation-gene stats table, keeping unmatched validation genes with null stats."""
    output_columns = [
        "val_gene",
        "ext_gene",
        "qval_diffexp",
        "pval_diffexp",
        "diffexp",
        "diffexp_se",
        "qval_dmr",
        "pval_dmr",
        "mean_methylation_difference",
        "germ_layer",
        "annotation_type",
    ]
    combined_df.drop_nulls("val_gene").select(output_columns).write_csv(
        snakemake.output["val_genes"], separator="\t"
    )


def add_val_genes(combined_df: pl.DataFrame) -> pl.DataFrame:
    """
    For each ext_gene in combined_df, look up whether it matches a known
    validation gene (or one of its synonyms), and add that match as a new
    `val_gene` column. No rows are removed or deduplicated; if ext_gene
    matches multiple validation genes, combined_df gains one row per match.
    """
    val_genes = pl.read_csv(
        snakemake.input["val_genes"], separator="\t", null_values="NA"
    )

    # Explode "LFS1,P53,TP53" into three rows (one synonym per row), so each
    # synonym can be matched against ext_gene with a plain equality join
    val_gene_synonyms = (
        val_genes.with_columns(pl.col("synonyms").str.split(","))
        .explode("synonyms")
        .rename({"synonyms": "synonym"})
    )

    combined_df = combined_df.join(
        val_gene_synonyms,
        left_on="ext_gene",
        right_on="synonym",
        how="left",
    )
    write_val_genes_table(combined_df)

    return combined_df


def plot_df(
    df: pl.DataFrame,
    x_domain: list[float],
    y_domain: list[float],
    title: str,
    show_axes: dict[str, bool],
) -> alt.Chart:
    """Scatter plot of diffexp vs. methylation difference, with a regression line and validation gene labels."""
    layer_select = alt.selection_point(
        fields=["germ_layer"], bind="legend", name="Germ Layer"
    )
    qval_slider = alt.param(
        name="qval_min",
        value=0.05,
        bind=alt.binding_range(min=0, max=1, step=0.01, name="q-value: "),
    )
    color = alt.condition(
        layer_select,
        alt.Color(
            "germ_layer:N",
            title="Germ Layer",
            scale=alt.Scale(
                domain=list(LAYER_COLORS.keys()), range=list(LAYER_COLORS.values())
            ),
        ),
        alt.value("lightgray"),
    )
    tooltip_cols = [
        "ext_gene",
        "diffexp",
        "mean_methylation_difference",
        "qval_combined",
        "qval_dmr",
        "qval_diffexp",
    ]

    base = alt.Chart(df.to_pandas(), title=title).transform_filter(
        alt.datum.qval_combined <= qval_slider
    )

    points = (
        base.mark_point(filled=True)
        .add_params(layer_select, qval_slider)
        .encode(
            x=alt.X(
                "diffexp:Q",
                scale=alt.Scale(domain=x_domain),
                title="Differential expression value" if show_axes.get("x") else None,
            ),
            y=alt.Y(
                "mean_methylation_difference:Q",
                scale=alt.Scale(domain=y_domain),
                title="DMR value" if show_axes.get("y") else None,
            ),
            size=alt.Size(
                "qval_combined:Q", title="q-value", scale=alt.Scale(range=[50, 1])
            ),
            opacity=alt.Opacity(
                "qval_combined:Q", scale=alt.Scale(range=[0.8, 0]), title="q-value"
            ),
            color=color,
            tooltip=tooltip_cols,
        )
    )

    regression_line = (
        base.transform_regression("diffexp", "mean_methylation_difference")
        .mark_line(size=2, color="blue")
        .encode(x="diffexp:Q", y="mean_methylation_difference:Q")
    )

    val_gene_labels = (
        base.transform_filter(alt.datum.val_gene != None)
        .transform_filter(layer_select)
        .mark_text(
            align="left",
            dx=3,
            dy=-3,
            fontWeight="bold",
            fontSize=10,
            color="black",
            stroke="white",
            strokeWidth=1,
        )
        .encode(
            x="diffexp:Q",
            y="mean_methylation_difference:Q",
            text="val_gene:N",
        )
    )

    return points + regression_line + val_gene_labels


diffexp_df = pl.concat(
    read_diffexp(path, layer, sign)
    for path, layer, sign in zip(diffexp_inputs, non_base_layers, diffexp_base_signs)
)
dmrs_df = pl.concat(
    read_dmrs(path, layer) for path, layer in zip(layer_inputs, non_base_layers)
)

combined_df = (
    diffexp_df.join(dmrs_df, on=["ext_gene", "germ_layer"], how="outer")
    .with_columns(
        pl.coalesce("ext_gene", "ext_gene_right").alias("ext_gene"),
        pl.coalesce("germ_layer", "germ_layer_right").alias("germ_layer"),
    )
    .drop("ext_gene_right", "germ_layer_right")
    .with_columns(
        qval_combined=pl.min_horizontal("qval_diffexp", "qval_dmr"),
        pval_combined=pl.min_horizontal("pval_diffexp", "pval_dmr"),
    )
)


combined_df = add_val_genes(combined_df)
annotation_type = None
if annotation_type != "unfiltered":
    combined_df = combined_df.filter(
        pl.col("annotation_type") == annotation_type
    ).drop_nulls(subset=pl.exclude("val_gene"))

# We want to have the same scale for all plots so we compute the domain here and not in the plotting function.
x_domain = [combined_df["diffexp"].min(), combined_df["diffexp"].max()]
y_domain = [
    combined_df["mean_methylation_difference"].min(),
    combined_df["mean_methylation_difference"].max(),
]

charts = [plot_df(combined_df, x_domain, y_domain, "All Layers", {"y": True})]
for layer in sorted(combined_df["germ_layer"].unique()):
    layer_df = combined_df.filter(pl.col("germ_layer") == layer)
    show_axes = {"x": layer in ("endoderm", "mesoderm"), "y": layer == "endoderm"}
    charts.append(plot_df(layer_df, x_domain, y_domain, layer, show_axes))

alt.concat(*charts, columns=2).save(snakemake.output.dmr_diffexp)


# This is only for downstream pathway analysis and is not used right now.
# diffexp_min, diffexp_max = symmetric_domain(combined_df["diffexp"])
# meth_diff_min, meth_diff_max = symmetric_domain(
#     combined_df["mean_methylation_difference"]
# )
# meth_diff_max_scaled = meth_diff_max * diffexp_max

# ranked_df = (
#     add_corner_distance_rank(combined_df, diffexp_max, meth_diff_max_scaled)
#     .select(
#         "ext_gene",
#         "ens_gene",
#         "germ_layer",
#         "qval_dmr",
#         "pval_dmr",
#         "diffexp",
#         "diffexp_se",
#         "qval_diffexp",
#         "pval_diffexp",
#         "qval_combined",
#         "pval_combined",
#         "mean_methylation_difference",
#         "ranked_meth_diffexp",
#     )
#     .with_row_index("row_id")
# )

combined_df.select(
    "ext_gene",
    "ens_gene",
    "germ_layer",
    "qval_dmr",
    "pval_dmr",
    "diffexp",
    "diffexp_se",
    "qval_diffexp",
    "pval_diffexp",
    "qval_combined",
    "pval_combined",
    "mean_methylation_difference",
    # "ranked_meth_diffexp",
).with_row_index("row_id").write_csv(snakemake.output.tsv, separator="\t")


# def symmetric_domain(series: pl.Series) -> tuple[float, float]:
#     """Return (-m, m) where m is the largest absolute value in `series`."""
#     bound = max(abs(series.min()), abs(series.max()))
#     return -bound, bound


# def add_corner_distance_rank(
#     df: pl.DataFrame, x_max: float, y_max: float
# ) -> pl.DataFrame:
#     """
#     Rank genes by how "extreme" their combined diffexp/methylation result is,
#     i.e. how close they sit to one of the plot's two outer corners
#     (top-left = down in expression & up in methylation, or the reverse).

#     Both axes are first scaled to comparable ranges, then for each point we
#     take the smaller of its distance to the top-left and bottom-right
#     corners. `ranked_meth_diffexp` is this distance converted into a single
#     signed score: positive and large for points near the top-left corner,
#     negative and large (more negative = further) for points near the
#     bottom-right corner. Sorting by |ranked_meth_diffexp| descending then
#     surfaces the most extreme genes in either direction first.
#     """
#     max_dist = (x_max**2 + y_max**2) ** 0.5

#     return (
#         df.with_columns(
#             mean_methylation_difference_scaled=pl.col("mean_methylation_difference")
#             * x_max
#         )
#         .with_columns(
#             dist_top_left=(
#                 (pl.col("diffexp") - (-x_max)) ** 2
#                 + (pl.col("mean_methylation_difference_scaled") - y_max) ** 2
#             ).sqrt(),
#             dist_bottom_right=(
#                 (pl.col("diffexp") - x_max) ** 2
#                 + (pl.col("mean_methylation_difference_scaled") - (-y_max)) ** 2
#             ).sqrt(),
#         )
#         .with_columns(
#             ranked_meth_diffexp=pl.when(
#                 pl.col("dist_top_left") < pl.col("dist_bottom_right")
#             )
#             .then(max_dist - pl.col("dist_top_left"))
#             .otherwise(-max_dist + pl.col("dist_bottom_right"))
#         )
#         .sort(pl.col("ranked_meth_diffexp").abs(), descending=True)
#     )
