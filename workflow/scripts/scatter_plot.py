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

    df[f"prob_present_{x_axis_name}_thresh"] = np.where(
        df[f"{x_axis_name}_prob_present"] >= 0.95, 0, 1
    )
    df[f"prob_present_{y_axis_name}_thresh"] = np.where(
        df[f"{y_axis_name}_prob_present"] >= 0.95, 0, 1
    )

    df["prob_present"] = (
        df[f"prob_present_{x_axis_name}_thresh"]
        + 2 * df[f"prob_present_{y_axis_name}_thresh"]
    )

    # Reassign category 4 for small methylation difference
    diff = abs(df[f"{x_axis_name}_methylation"] - df[f"{y_axis_name}_methylation"])
    df["prob_present"] = np.where(diff <= 20, 4, df["prob_present"])

    color_map = {
        0: "black",
        1: "red",
        2: "blue",
        3: "purple",
        4: "grey",
    }

    category_desc = {
        0: "No prob present",
        1: f"Only {y_axis_name} significant",
        2: f"Only {x_axis_name} significant",
        3: "Both significant",
        4: "Beta val <= 20",
    }

    df["color"] = df["prob_present"].map(color_map)
    df["category"] = df["prob_present"].map(category_desc)
    print(df)
    chart = (
        alt.Chart(df)
        .mark_circle(size=30, opacity=0.3)
        .encode(
            x=alt.X(f"{x_axis_name}_methylation", title=f"{x_axis_name} Methylation"),
            y=alt.Y(f"{y_axis_name}_methylation", title=f"{y_axis_name} Methylation"),
            color=alt.Color(
                "category:N",
                title="Prob Present Categories",
                scale=alt.Scale(
                    domain=list(category_desc.values()), range=list(color_map.values())
                ),
            ),
            tooltip=[
                f"{x_axis_name}_methylation",
                f"{y_axis_name}_methylation",
                "category",
            ],
        )
        .properties(
            title=f"{x_axis_name} vs. {y_axis_name} — RMSE: {rmse:.2f}, N: {len(df)}",
            width=500,
            height=500,
        )
        .configure_legend(titleFontSize=12, labelFontSize=10)
    )

    chart.save(output)


# Einstellungen für DataFrames
alt.data_transformers.enable("vegafusion")
pd.set_option("display.max_columns", None)

df = pd.read_parquet(snakemake.input["calls"], engine="pyarrow")
x_axis = snakemake.params["group1"]
y_axis = snakemake.params["group2"]
df = df[(df[f"{x_axis}_coverage"] > 50) & (df["{y_axis}_coverage"] > 50)]


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
