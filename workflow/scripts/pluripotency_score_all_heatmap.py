import sys

import altair as alt
import numpy as np
import polars as pl

# Maps the methylation column names in the input dataframe to short layer names.
METH_TO_LAYER = {
    "psc_methylation": "psc",
    "endoderm_methylation": "endoderm",
    "mesoderm_methylation": "mesoderm",
    "ectoderm_methylation": "ectoderm",
}

# For each biomarker, which germ layer(s) it is expected to mark ("target" vs "other").
# "endomeso" markers are shared between endoderm and mesoderm.
BIOMARKER_TO_TARGET = {
    "psc": {"psc"},
    "endoderm": {"endoderm"},
    "mesoderm": {"mesoderm"},
    "ectoderm": {"ectoderm"},
    "endomeso": {"endoderm", "mesoderm"},
}

LAYER_COLORS = {
    "psc": "#05AA8F",
    "endoderm": "#D81B60",
    "mesoderm": "#1E88E5",
    "ectoderm": "#FFC107",
}

BIOMARKER_POSITIONS = pl.DataFrame(
    [
        # PSC
        {
            "chromosome": "2",
            "position": 15938891,
            "direction": -1,
            "cg_id": "cg21699252",
            "biomarker": "psc",
        },
        {
            "chromosome": "10",
            "position": 35594676,
            "direction": -1,
            "cg_id": "cg00933813",
            "biomarker": "psc",
        },
        {
            "chromosome": "4",
            "position": 168841140,
            "direction": -1,
            "cg_id": "cg00661673",
            "biomarker": "psc",
        },
        # Endoderm
        {
            "chromosome": "6",
            "position": 12886978,
            "direction": 1,
            "cg_id": "cg20548013",
            "biomarker": "endoderm",
        },
        {
            "chromosome": "11",
            "position": 8840472,
            "direction": 1,
            "cg_id": "cg14521421",
            "biomarker": "endoderm",
        },
        {
            "chromosome": "8",
            "position": 125637563,
            "direction": -1,
            "cg_id": "cg08913523",
            "biomarker": "endoderm",
        },
        # Mesoderm
        {
            "chromosome": "2",
            "position": 128638846,
            "direction": 1,
            "cg_id": "cg14708360",
            "biomarker": "mesoderm",
        },
        {
            "chromosome": "17",
            "position": 15966293,
            "direction": 1,
            "cg_id": "cg08826152",
            "biomarker": "mesoderm",
        },
        {
            "chromosome": "12",
            "position": 122872581,
            "direction": 1,
            "cg_id": "cg11599718",
            "biomarker": "mesoderm",
        },
        # Ectoderm
        {
            "chromosome": "15",
            "position": 71314056,
            "direction": -1,
            "cg_id": "cg01907071",
            "biomarker": "ectoderm",
        },
        {
            "chromosome": "5",
            "position": 107661548,
            "direction": -1,
            "cg_id": "cg18118164",
            "biomarker": "ectoderm",
        },
        {
            "chromosome": "14",
            "position": 68569167,
            "direction": -1,
            "cg_id": "cg13075942",
            "biomarker": "ectoderm",
        },
        # Endomesoderm
        {
            "chromosome": "5",
            "position": 111309143,
            "direction": 1,
            "cg_id": "cg23385847",
            "biomarker": "endomeso",
        },
        {
            "chromosome": "5",
            "position": 24208809,
            "direction": 1,
            "cg_id": "cg24919344",
            "biomarker": "endomeso",
        },
        {
            "chromosome": "10",
            "position": 33773306,
            "direction": 1,
            "cg_id": "cg11147278",
            "biomarker": "endomeso",
        },
    ],
    schema={
        "chromosome": pl.Utf8,
        "position": pl.Int64,
        "direction": pl.Int64,
        "cg_id": pl.Utf8,
        "biomarker": pl.Utf8,
    },
)


def adjust_for_direction_expr(methylation_col, direction_col):
    """
    When the position is defined as hypomethylated we return the complement of the methylation value.
    """
    return (
        pl.when(direction_col == 1).then(methylation_col).otherwise(1 - methylation_col)
    )


def filter_on_biomarker(df, biomarker_positions):
    """
    Filter the input data down to the known biomarker CpG positions
    """
    matched = df.join(
        biomarker_positions,
        on=["chromosome", "position"],
        how="inner",
    )
    long_df = matched.unpivot(
        index=["chromosome", "position", "direction", "biomarker", "cg_id"],
        on=list(METH_TO_LAYER.keys()),
        variable_name="layer",
        value_name="methylation",
    )
    long_df = long_df.with_columns(
        pl.col("layer").replace(METH_TO_LAYER),
        (pl.col("methylation") / 100),
    )
    long_df = long_df.with_columns(
        adjust_for_direction_expr(pl.col("methylation"), pl.col("direction")).alias(
            "adjusted_methylation"
        )
    )
    long_df = long_df.with_columns(
        pl.when(
            pl.struct(["biomarker", "layer"]).map_elements(
                lambda r: r["layer"] in BIOMARKER_TO_TARGET[r["biomarker"]],
                return_dtype=pl.Boolean,
            )
        )
        .then(pl.lit("target"))
        .otherwise(pl.lit("other"))
        .alias("type")
    )

    return long_df


def build_biomarker_subtitles(long_df):
    """
    Build a per-biomarker subtitle string summarizing the summed adjusted
    methylation per germ layer
    """
    score_df = (
        long_df.group_by(["biomarker", "layer"])
        .agg(pl.col("adjusted_methylation").sum())
        .sort(["biomarker", "layer"])
    )

    subtitles = {}
    for biomarker, group in score_df.group_by("biomarker"):
        # group_by on a single key still returns a 1-tuple key
        biomarker = biomarker[0]
        parts = [
            f"{row['layer'].replace('derm', '')}: {row['adjusted_methylation']:.2f}"
            for row in group.iter_rows(named=True)
        ]
        subtitles[biomarker] = ", ".join(parts)
    return subtitles


def make_biomarker_chart(df_sub, biomarker, subtitle, show_y_axis, is_last_column):
    """Build a single point chart (methylation vs. cell-fate type) for one biomarker."""
    # Add jitter
    rng = np.random.default_rng(42)
    jitter = rng.normal(loc=0.0, scale=0.06, size=df_sub.height)
    df_sub = df_sub.with_columns(pl.Series("y_jitter", jitter))
    return (
        alt.Chart(df_sub.to_pandas())
        .mark_point(size=100, filled=True)
        .encode(
            x=alt.X(
                "adjusted_methylation:Q",
                scale=alt.Scale(domain=[0, 1]),
                title="methylation" if is_last_column else None,
                axis=alt.Axis(labels=True, ticks=True, domain=True),
            ),
            y=alt.Y(
                "type:N",
                title="cell fate" if show_y_axis else None,
                axis=alt.Axis(labels=show_y_axis, ticks=False, domain=True),
            ),
            yOffset=alt.YOffset(
                "y_jitter:Q",
                scale=alt.Scale(domain=[-1, 1]),
            ),
            color=alt.Color(
                "layer:N",
                title="Germ Layer",
                scale=alt.Scale(
                    domain=list(LAYER_COLORS.keys()),
                    range=list(LAYER_COLORS.values()),
                ),
            ),
        )
        .properties(
            title={
                "text": biomarker,
                "subtitle": subtitle,
                "subtitleFontSize": 9,
                "subtitleFontStyle": "Liberation Sans",
            },
        )
    )


sys.stderr = open(snakemake.log[0], "w", buffering=1)
pl.Config.set_tbl_cols(
    -1
)  # show all columns when printing, like pandas display.max_columns
df = pl.read_parquet(snakemake.input[0]).with_columns(
    pl.col("chromosome").cast(pl.Utf8)
)

df = filter_on_biomarker(df, BIOMARKER_POSITIONS)
subtitles = build_biomarker_subtitles(df)

charts = []
biomarkers = ["psc"] + df["biomarker"].unique(maintain_order=True).to_list()
ncols = 2
for i, biomarker in enumerate(biomarkers):
    df_sub = df.filter(pl.col("biomarker") == biomarker)
    show_y_axis = i % ncols == 0
    is_last_column = i == len(biomarkers) - 1
    print(df_sub)
    charts.append(
        make_biomarker_chart(
            df_sub,
            biomarker,
            subtitles.get(biomarker, ""),
            show_y_axis,
            is_last_column,
        )
    )

chart = alt.concat(*charts, columns=ncols)
chart.save(snakemake.output[0])
