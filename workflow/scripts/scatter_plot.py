import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import altair as alt
import os

sys.stderr = open(snakemake.log[0], "w", buffering=1)


def compute_rmse(df, x_axis_name, y_axis_name):
    squared_errors = (
        df[f"{x_axis_name}_methylation"] - df[f"{y_axis_name}_methylation"]
    ) ** 2
    mean_squared_error = squared_errors.mean()
    return np.sqrt(mean_squared_error)


def plot_meth_vals(df, output, x_axis_name, y_axis_name):
    rmse = compute_rmse(df, x_axis_name, y_axis_name)

    if snakemake.params["meth_caller"] == "varlo":
        df[f"prob_present_{x_axis_name}_thresh"] = np.where(
            df[f"{x_axis_name}_prob_present"] >= 0.95, 1, 0
        )
        df[f"prob_present_{y_axis_name}_thresh"] = np.where(
            df[f"{y_axis_name}_prob_present"] >= 0.95, 1, 0
        )
        # 0 if both are significant, 1 if only x is significant, 2 if only y is significant, 3 if neither
        df["prob_present"] = (
            df[f"prob_present_{x_axis_name}_thresh"]
            + 2 * df[f"prob_present_{y_axis_name}_thresh"]
        )
        df = df[df["prob_present"] != 0]

        diff = abs(df[f"{x_axis_name}_methylation"] - df[f"{y_axis_name}_methylation"])
        df["prob_present"] = np.where(diff <= 20, 0, df["prob_present"])

        color_map = {
            0: "grey",
            1: "red",
            2: "blue",
            3: "black",
        }

        category_desc = {
            0: "Beta val <= 20",
            1: f"Only {x_axis_name} significant",
            2: f"Only {y_axis_name} significant",
            3: "Both significant",
        }

        df["color"] = df["prob_present"].map(color_map)
        df["category"] = df["prob_present"].map(category_desc)

    else:
        df["color"] = "grey"
        df["category"] = "No probs given"

    print(df, file=sys.stderr)
    chart = (
        alt.Chart(df)
        .mark_circle(size=15, opacity=0.5)
        .encode(
            x=alt.X(f"{x_axis_name}_methylation", title=f"{x_axis_name} Methylation"),
            y=alt.Y(f"{y_axis_name}_methylation", title=f"{y_axis_name} Methylation"),
            color=alt.Color(
                "category:N",
                title="Prob Present Categories",
                scale=alt.Scale(
                    domain=df["category"].unique().tolist(),
                    range=df["color"].unique().tolist(),
                ),
            ),
            tooltip=[
                "chromosome",
                "position",
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


alt.data_transformers.enable("vegafusion")
pd.set_option("display.max_columns", None)


df = pd.read_parquet(snakemake.input["calls"], engine="pyarrow")

print(df, file=sys.stderr)
x_axis = snakemake.params["group2"]
y_axis = snakemake.params["group1"]

df = df[(df[f"{x_axis}_coverage"] > 50) & (df[f"{y_axis}_coverage"] > 50)]
print(df, file=sys.stderr)


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
