import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import altair as alt
import os


sys.stderr = open(snakemake.log[0], "w", buffering=1)


def compute_rmse(df, cell_type):
    squared_errors = (
        df[f"{cell_type}_methylation_np"] - df[f"{cell_type}_methylation_pb"]
    ) ** 2
    mean_squared_error = squared_errors.mean()
    return np.sqrt(mean_squared_error)


def plot_meth_vals(df, output, cell_type):
    rmse = compute_rmse(df, cell_type)
    df[f"prob_present_{cell_type}_thresh_np"] = np.where(
        df[f"{cell_type}_prob_present_np"] >= 0.95, 1, 0
    )
    df[f"prob_present_{cell_type}_thresh_pb"] = np.where(
        df[f"{cell_type}_prob_present_pb"] >= 0.95, 1, 0
    )
    df["prob_present"] = (
        df[f"prob_present_{cell_type}_thresh_np"]
        + 2 * df[f"prob_present_{cell_type}_thresh_pb"]
    )

    df = df[df["prob_present"] != 0]


    x = df[f"{cell_type}_methylation_np"]
    y = df[f"{cell_type}_methylation_pb"]

    color_map = {4: "grey", 3: "purple", 2: "blue", 1: "red", 0: "black"}
    category_desc = {
        0: "no prob present",
        1: f"Only {cell_type} in pb significant",
        2: f"Only {cell_type} in np significant",
        3: "Both significant",
        4: "Beta val <= 20",
    }


    color_map = {
        0: "grey",
        1: "red",
        2: "blue",
        3: "black",
    }

    category_desc = {
        0: "Beta val <= 20",
        1: f"Only {cell_type} in np significant",
        2: f"Only {cell_type} in pb significant",
        3: "Both significant",
    }


    df["prob_present"] = np.where(
        abs(df[f"{cell_type}_methylation_np"] - df[f"{cell_type}_methylation_pb"])
        <= 20,
        0,
        df["prob_present"],
    )
    df["color"] = df["prob_present"].map(color_map)


    df["category"] = df["prob_present"].map(category_desc)
    print(df)
    chart = (
        alt.Chart(df)
        .mark_circle(size=15, opacity=0.5)
        .encode(
            x=alt.X(f"{cell_type}_methylation_np", title=f"{cell_type} Nanopore Methylation"),
            y=alt.Y(f"{cell_type}_methylation_pb", title=f"{cell_type} PacBio Methylation"),
            color=alt.Color(
                "category:N",
                title="Prob Present Categories",
                scale=alt.Scale(
                    domain=list(category_desc.values()), range=list(color_map.values())
                ),
            ),
        )
        .properties(
            title=f"{cell_type} Nanopore vs. {cell_type} PacBio\nRmse: {rmse:.2f}, Datapoints: {len(df)}",
            width=500,
            height=500,
        )
        .configure_legend(titleFontSize=12, labelFontSize=10)
    )

    chart.save(output)


alt.data_transformers.enable("vegafusion")
pd.set_option("display.max_columns", None)

df = pd.read_parquet(snakemake.input["pacbio"], engine="pyarrow")
df_pb = df[
    (df["mesoderm_coverage"] > 50)
    & (df["endoderm_coverage"] > 50)
]
df = pd.read_parquet(snakemake.input["nanopore"], engine="pyarrow")
df_np = df[
    (df["mesoderm_coverage"] > 50)
    & (df["endoderm_coverage"] > 50)
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
