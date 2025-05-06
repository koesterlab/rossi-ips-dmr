import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import altair as alt
import os


def compute_rmse(df, x_axis_name, y_axis_name):
    print(df, x_axis_name)
    squared_errors = (
        df[f"{x_axis_name}_methylation"] - df[f"{y_axis_name}_methylation"]
    ) ** 2
    mean_squared_error = squared_errors.mean()
    return np.sqrt(mean_squared_error)


def plot_meth_vals(df, output, x_axis_name, y_axis_name):
    rmse = compute_rmse(df, x_axis_name, y_axis_name)
    df[f"prob_present_{x_axis_name}_thresh"] = np.where(df[f"{x_axis_name}_prob_present"] >= 0.95, 0, 1)
    df[f"prob_present_{y_axis_name}_thresh"] = np.where(df[f"{y_axis_name}_prob_present"] >= 0.95, 0, 1)
    df["prob_present"] = df[f"prob_present_{x_axis_name}_thresh"] + 2 * df[f"prob_present_{y_axis_name}_thresh"]

    x = df[f"{x_axis_name}_methylation"]
    y = df[f"{y_axis_name}_methylation"]

    color_map = {4: "grey", 3: "purple", 2: "blue", 1: "red", 0: "black"}
    df["color"] = df["prob_present"].map(color_map)
    category_desc = {
        0: "no prob present",
        1: f"Only {y_axis_name} significant",
        2: f"Only {x_axis_name} significant",
        3: "Both significant",
        4: "Beta val <= 20",
    }

    df["prob_present"] = np.where(
        abs(df[f"{x_axis_name}_methylation"] - df[f"{y_axis_name}_methylation"]) <= 20,
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

    plt.xlabel(f"{x_axis_name} Methylation")
    plt.ylabel(f"{y_axis_name} Methylation")
    plt.title(
        f"{x_axis_name} vs. {y_axis_name}\nRmse: {rmse:.2f}, Datapoints: {len(df)}"
    )

    plt.legend(title="Prob Present Categories")

    plt.tight_layout()
    plt.savefig(output, dpi=200)
    plt.close()


# Einstellungen für DataFrames
alt.data_transformers.enable("vegafusion")
pd.set_option("display.max_columns", None)

df = pd.read_parquet(snakemake.input["calls"], engine="pyarrow")
df = df[(df['psc_coverage'] > 50) & (df['mesoderm_coverage'] > 50) & (df['endoderm_coverage'] > 50) & (df['ectoderm_coverage'] > 50)]
x_axis = snakemake.params["group1"]
y_axis = snakemake.params["group2"]


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


plot_meth_vals(df, snakemake.output[0], x_axis, y_axis)
