import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import altair as alt
import os


def compute_rmse(df, cell_type):
    squared_errors = (
        df[f"{cell_type}_methylation_np"] - df[f"{cell_type}_methylation_pb"]
    ) ** 2
    mean_squared_error = squared_errors.mean()
    return np.sqrt(mean_squared_error)


def plot_meth_vals(df, output, cell_type):
    rmse = compute_rmse(df, cell_type)
    df[f"prob_present_{cell_type}_thresh_np"] = np.where(
        df[f"{cell_type}_prob_present_np"] >= 0.95, 0, 1
    )
    df[f"prob_present_{cell_type}_thresh_pb"] = np.where(
        df[f"{cell_type}_prob_present_pb"] >= 0.95, 0, 1
    )
    df["prob_present"] = (
        df[f"prob_present_{cell_type}_thresh_np"]
        + 2 * df[f"prob_present_{cell_type}_thresh_pb"]
    )

    x = df[f"{cell_type}_methylation_np"]
    y = df[f"{cell_type}_methylation_pb"]

    color_map = {4: "grey", 3: "purple", 2: "blue", 1: "red", 0: "black"}
    df["color"] = df["prob_present"].map(color_map)
    category_desc = {
        0: "no prob present",
        1: f"Only {cell_type} in pb significant",
        2: f"Only {cell_type} in np significant",
        3: "Both significant",
        4: "Beta val <= 20",
    }

    df["prob_present"] = np.where(
        abs(df[f"{cell_type}_methylation_np"] - df[f"{cell_type}_methylation_pb"])
        <= 20,
        4,
        df["prob_present"],
    )

    for category, color in color_map.items():
        subset = df[df["prob_present"] == category]
        plt.scatter(
            subset[x.name],
            subset[y.name],
            color=color,
            label=f"{category_desc[category]}",
            s=3,
            alpha=0.3,
        )

    plt.xlabel(f"{cell_type} Nanopore Methylation")
    plt.ylabel(f"{cell_type} PacBio Methylation")
    plt.title(
        f"{cell_type} Nanopore vs. {cell_type} PacBio\nRmse: {rmse:.2f}, Datapoints: {len(df)}"
    )

    plt.legend(title="Prob Present Categories")

    plt.tight_layout()
    plt.savefig(output, dpi=200)
    plt.close()


# Einstellungen für DataFrames
alt.data_transformers.enable("vegafusion")
pd.set_option("display.max_columns", None)

df = pd.read_parquet(snakemake.input["pacbio"], engine="pyarrow")
df_pb = df[
    (df["psc_coverage"] > 50)
    & (df["mesoderm_coverage"] > 50)
    & (df["endoderm_coverage"] > 50)
    & (df["ectoderm_coverage"] > 50)
]
df = pd.read_parquet(snakemake.input["nanopore"], engine="pyarrow")
df_np = df[
    (df["psc_coverage"] > 50)
    & (df["mesoderm_coverage"] > 50)
    & (df["endoderm_coverage"] > 50)
    & (df["ectoderm_coverage"] > 50)
]
cell_type = snakemake.params["group"]
df = df_np.merge(df_pb, on=["chromosome", "position"], suffixes=("_np", "_pb"))
# y_axis = snakemake.params["group2"]


# most_variable_positions_df = pd.read_parquet(
#     snakemake.input["most_variable_positions"], engine="pyarrow"
# )
# differentiated_filtered = pd.merge(
#     differentiated_df,
#     most_variable_positions_df[["chromosome", "position"]],
#     on=["chromosome", "position"],
#     how="inner",
# )
# undifferentiated_filtered = pd.merge(
#     undifferentiated_df,
#     most_variable_positions_df[["chromosome", "position"]],
#     on=["chromosome", "position"],
#     how="inner",
# )


plot_meth_vals(df, snakemake.output[0], cell_type)
