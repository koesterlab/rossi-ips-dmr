import sys

import altair as alt
import numpy as np
import pandas as pd

sys.stderr = open(snakemake.log[0], "w", buffering=1)
alt.data_transformers.enable("vegafusion")

MAX_POINTS = 100_000
# Minimum posterior probability for a call to count as significant
PROB_THRESHOLD = 0.95
# Absolute methylation differences (in %) up to this value are shown in grey
MIN_METH_DIFF = 20

CATEGORY_COLORS = {0: "grey", 1: "red", 2: "blue", 3: "black"}


def compute_rmse(df, x_col, y_col):
    return np.sqrt(((df[x_col] - df[y_col]) ** 2).mean())


def categorize(df, x_axis, y_axis):
    """
    Encode which call is significant: 1 = only x, 2 = only y, 3 = both.
    Sites where neither call is significant are removed and sites with a small
    methylation difference are put into category 0.
    """
    category = (df[f"{x_axis}_prob_identified"] >= PROB_THRESHOLD).astype(int) + 2 * (
        df[f"{y_axis}_prob_identified"] >= PROB_THRESHOLD
    ).astype(int)
    df = df[category != 0].copy()
    small_diff = (df[f"{x_axis}_methylation"] - df[f"{y_axis}_methylation"]).abs() <= MIN_METH_DIFF
    df["prob_present"] = np.where(small_diff, 0, category[category != 0])
    descriptions = {
        0: f"Beta val <= {MIN_METH_DIFF}",
        1: f"Only {x_axis} significant",
        2: f"Only {y_axis} significant",
        3: "Both significant",
    }
    df["category"] = df["prob_present"].map(descriptions)
    df["color"] = df["prob_present"].map(CATEGORY_COLORS)
    return df


def plot_meth_vals(df, output, x_axis, y_axis):
    x_col, y_col = f"{x_axis}_methylation", f"{y_axis}_methylation"
    rmse = compute_rmse(df, x_col, y_col)
    df = categorize(df, x_axis, y_axis)
    # Domain and range must be in the same order
    legend = df[["category", "color"]].drop_duplicates()

    chart = (
        alt.Chart(df)
        .mark_circle(size=15, opacity=0.5)
        .encode(
            x=alt.X(x_col, title=f"{x_axis} Methylation"),
            y=alt.Y(y_col, title=f"{y_axis} Methylation"),
            color=alt.Color(
                "category:N",
                title="Prob Present Categories",
                scale=alt.Scale(
                    domain=legend["category"].tolist(),
                    range=legend["color"].tolist(),
                ),
            ),
        )
        .properties(
            title=f"{x_axis} vs. {y_axis} — RMSE: {rmse:.2f}, N: {len(df)}",
            width=500,
            height=500,
        )
        .configure_legend(titleFontSize=12, labelFontSize=10)
    )
    chart.save(output)


x_axis = snakemake.params["group2"]
y_axis = snakemake.params["group1"]

df = pd.read_parquet(
    snakemake.input["calls"],
    engine="pyarrow",
    columns=[
        f"{x_axis}_methylation",
        f"{y_axis}_methylation",
        f"{x_axis}_prob_identified",
        f"{y_axis}_prob_identified",
    ],
)
df = df.sample(min(len(df), MAX_POINTS), random_state=42)
plot_meth_vals(df, snakemake.output[0], x_axis, y_axis)
