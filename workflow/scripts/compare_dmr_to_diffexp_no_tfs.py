import sys
from unittest.main import main

import altair as alt
import polars as pl

sys.stderr = open(snakemake.log[0], "w", buffering=1)

pl.Config.set_tbl_rows(10000000)
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
non_base_layers: list[str] = snakemake.params["non_base_layers"]


diffexp_inputs = snakemake.input.diffexp
diffexp_base_signs = snakemake.params["diffexp_base_signs"]
layer_inputs = [snakemake.input.layer1, snakemake.input.layer2, snakemake.input.layer3]


def read_dmrs(path: str, layer: str) -> pl.DataFrame:
    df = (
        pl.read_csv(
            path,
            separator="\t",
            schema_overrides={"chr": pl.Utf8},
            null_values="NA",
        )
        .with_columns(pl.lit(layer).alias("germ_layer"))
        .filter(pl.col("annotation_type") == annotation_type)
        .select(
            "ext_gene",
            "transcriptId",
            "ens_gene",
            "annotation_type",
            "mean_methylation_difference",
            "absolute_signed_pi_val",
            "qval",
            "pval",
            "chr",
            "start_dmr",
            "end_dmr",
            # "annotation",
            "germ_layer",
        )
        .rename(
            {
                "qval": "qval_dmr",
                "pval": "pval_dmr",
            }
        )
    )
    print(df, file=sys.stderr)

    return df


def extract_comparison_name(filename: str) -> str:
    """
    Extract the comparison name from diffexp filename.
    e.g., 'ectoderm_vs_psc_new' from 'ectoderm_vs_psc_new.genes-representative.diffexp_postprocessed.tsv'
    """
    # Remove path and extension
    basename = filename.split("/")[-1]
    # Remove the '.genes-representative.diffexp_postprocessed.tsv' suffix
    comparison = basename.replace(".genes-representative.diffexp_postprocessed.tsv", "")
    return comparison


def read_diffexp(path: str, layer: str, sign: int) -> pl.DataFrame:
    """
    Read diffexp file and rename qval/pval columns to include the comparison name.
    Multiply beta columns by the sign parameter.
    """
    df = pl.read_csv(path, separator="\t", null_values="NA")

    # Rename qval and pval columns to include the comparison name

    # Multiply beta columns by sign
    effect_col = [
        col
        for col in df.columns
        if col.endswith("-") or col.endswith("+") and col.startswith("b_")
    ][0]
    signed_pi_col = [
        col
        for col in df.columns
        if col.endswith("-") or col.endswith("+") and col.startswith("signed_pi")
    ][0]
    se_col = [col for col in df.columns if col.endswith("_se")][0]
    df = (
        df.with_columns((pl.col(effect_col) * sign).alias("diffexp"))
        .with_columns((pl.col(signed_pi_col) * sign).alias("signed_pi_col"))
        .with_columns((pl.lit(layer).alias("germ_layer")))
        .with_columns((pl.col(se_col)).alias("diffexp_se"))
        .select(
            [
                "ext_gene",
                "ens_gene",
                "target_id",
                "qval",
                "pval",
                "diffexp",
                "diffexp_se",
                "signed_pi_col",
                "germ_layer",
            ]
        )
        .rename(
            {
                "qval": "qval_diffexp",
                "pval": "pval_diffexp",
            }
        )
    )

    return df


def find_val_genes(diffexp_df: pl.DataFrame, dmrs_df: pl.DataFrame) -> pl.DataFrame:
    val_genes = pl.read_csv(
        snakemake.input["val_genes"], separator="\t", null_values="NA"
    )
    print(diffexp_df.filter(pl.col("ext_gene") == "FEZF1"))
    print(dmrs_df.filter(pl.col("ext_gene") == "FEZF1"))
    common_df = diffexp_df.join(dmrs_df, on="ext_gene", how="outer")
    common_df = common_df.with_columns(
        pl.coalesce("ext_gene", "ext_gene_right").alias("ext_gene"),
        pl.coalesce("germ_layer", "germ_layer_right").alias("germ_layer"),
    ).drop("ext_gene_right", "germ_layer_right")
    print(common_df.filter(pl.col("ext_gene") == "FEZF1"))

    result = val_genes.join(common_df, on="ext_gene", how="left").select(
        "ext_gene",
        "qval_diffexp",
        "pval_diffexp",
        "diffexp",
        "diffexp_se",
        "qval_dmr",
        "pval_dmr",
        "mean_methylation_difference",
        # "absolute_signed_pi_val",
        "germ_layer",
        # "annotation_type",
    )

    # Write output
    #  WidBk?
    #
    result.write_csv(snakemake.output["val_genes"], separator="\t")


def merge_diffexp_dmr(diffexp_df: pl.DataFrame, dmrs_df: pl.DataFrame) -> pl.DataFrame:
    common_df = diffexp_df.join(
        dmrs_df, on=["ext_gene", "germ_layer"], how="inner"
    ).filter(pl.col("annotation_type") == annotation_type)

    common_df = common_df.with_columns(
        pl.min_horizontal(pl.col("qval_diffexp"), pl.col("qval_dmr")).alias(
            "qval_combined"
        )
    ).with_columns(
        pl.min_horizontal(pl.col("qval_diffexp"), pl.col("pval_dmr")).alias(
            "pval_combined"
        )
    )
    return common_df


def plot_df(common_df: pl.DataFrame):
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
        name="qval_min",
        value=0.05,
        bind=alt.binding_range(
            min=0, max=1, step=0.01, name="min(qval_dmr, qval_diffexp): "
        ),
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

    # Prepare tooltip columns - include all qval and pval columns that exist
    tooltip_cols = [
        "ext_gene",
        "diffexp",
        "mean_methylation_difference",
        "qval_combined",
        "qval_dmr",
        "qval_diffexp",
    ]

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
            tooltip=tooltip_cols,
            color=color,
            opacity=alt.Opacity(
                "qval_combined:Q",
                scale=alt.Scale(range=[0.5, 0]),
                title="max(qval1, qval2)",
            ),
        )
        .add_params(layer_select, qval_slider)
    )

    chart.save(snakemake.output.html)


if __name__ == "__main__":
    diffexp_df = pl.concat(
        [
            read_diffexp(path, layer, sign)
            for path, layer, sign in zip(
                diffexp_inputs, non_base_layers, diffexp_base_signs
            )
        ]
    )
    dmrs_df = pl.concat(
        [read_dmrs(path, layer) for path, layer in zip(layer_inputs, non_base_layers)]
    )

    find_val_genes(diffexp_df, dmrs_df)

    common_df = merge_diffexp_dmr(diffexp_df, dmrs_df)

    plot_df(common_df)

    diffexp_min, diffexp_max, meth_diff_min, meth_diff_max = common_df.select(
        pl.col("diffexp").min().alias("x_min"),
        pl.col("diffexp").max().alias("x_max"),
        pl.col("mean_methylation_difference").min().alias("y_min"),
        pl.col("mean_methylation_difference").max().alias("y_max"),
    ).row(0)
    # The values should be symmetric around zero
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
                    + (
                        pl.col("mean_methylation_difference_scaled")
                        - meth_diff_max_scaled
                    )
                    ** 2
                ).sqrt()
            ).alias("dist_top_left"),
            (
                (
                    (pl.col("diffexp") - diffexp_max) ** 2
                    + (
                        pl.col("mean_methylation_difference_scaled")
                        - meth_diff_min_scaled
                    )
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
                "qval_dmr",
                "pval_dmr",
                "diffexp",
                "diffexp_se",
                "qval_diffexp",
                "pval_diffexp",
                "qval_combined",
                "pval_combined",
                "mean_methylation_difference",
                "ranked_meth_diffexp",
            ]
        )
        .with_row_index("row_id")
    )

    common_df.write_csv(snakemake.output.tsv, separator="\t")
