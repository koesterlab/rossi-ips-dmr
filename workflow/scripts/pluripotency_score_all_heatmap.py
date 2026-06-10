import sys

import altair as alt
import pandas as pd


def calculate_adjusted_methylation(methylation_value, direction):
    """
    Adjust methylation based on direction:
    - If direction=1: return methylation as is
    - If direction=-1: return 1 - methylation
    """
    return methylation_value if direction == 1 else 1 - methylation_value


sys.stderr = open(snakemake.log[0], "w", buffering=1)

pd.set_option("display.max_columns", None)
df = pd.read_parquet(snakemake.input, engine="pyarrow")

METHYLATION_COLS = {
    "psc_methylation": "psc",
    "endoderm_methylation": "endoderm",
    "mesoderm_methylation": "mesoderm",
    "ectoderm_methylation": "ectoderm",
}

BIOMARKER_POSITIONS = pd.DataFrame(
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
    ]
)

BIOMARKER_POSITIONS["chromosome"] = BIOMARKER_POSITIONS["chromosome"].astype(str)

BIOMARKER_TO_TARGET = {
    "psc": {"psc"},
    "endoderm": {"endoderm"},
    "mesoderm": {"mesoderm"},
    "ectoderm": {"ectoderm"},
    "endomeso": {"endoderm", "mesoderm"},
}

LAYER_COLORS = {
    "psc": "#1f77b4",
    "endoderm": "#ff7f0e",
    "mesoderm": "#2ca02c",
    "ectoderm": "#d62728",
}

df["chromosome"] = df["chromosome"].astype(str)
df["position_pair"] = list(zip(df["chromosome"], df["position"]))

biomarker_pairs = set(
    zip(BIOMARKER_POSITIONS["chromosome"], BIOMARKER_POSITIONS["position"])
)

long_df = (
    df[df["position_pair"].isin(biomarker_pairs)]
    .copy()
    .merge(
        BIOMARKER_POSITIONS[
            ["chromosome", "position", "direction", "biomarker", "cg_id"]
        ],
        on=["chromosome", "position"],
        how="left",
    )
    .melt(
        id_vars=["chromosome", "position", "direction", "biomarker", "cg_id"],
        value_vars=list(METHYLATION_COLS.keys()),
        var_name="layer",
        value_name="methylation",
    )
)

long_df["layer"] = long_df["layer"].map(METHYLATION_COLS)
long_df["methylation"] = long_df["methylation"] / 100


def adjust(m, d):
    return m if d == 1 else 1 - m


long_df["adjusted_methylation"] = long_df.apply(
    lambda r: adjust(r["methylation"], r["direction"]), axis=1
)


long_df["type"] = long_df.apply(
    lambda r: (
        "target" if r["layer"] in BIOMARKER_TO_TARGET[r["biomarker"]] else "other"
    ),
    axis=1,
)


score_df = (
    long_df.groupby(["biomarker", "layer"])["adjusted_methylation"].sum().reset_index()
)

score_map = (
    score_df.groupby("biomarker")
    .apply(
        lambda x: ", ".join(
            f"{row.layer.replace('derm', '')}: {row.adjusted_methylation:.2f}"
            for row in x.itertuples()
        )
    )
    .to_dict()
)

long_df["subtitle"] = long_df["biomarker"].map(score_map)


charts = []
subtitles = {}
ncols = 2

for i, biomarker in enumerate(long_df["biomarker"].unique()):
    df_sub = long_df[long_df["biomarker"] == biomarker]
    print(df_sub)
    show_y = i % ncols == 0
    # show_x = i == 4
    show_x = True

    c = (
        alt.Chart(df_sub)
        .mark_point(size=100, filled=True)
        .encode(
            x=alt.X(
                "methylation:Q",
                scale=alt.Scale(domain=[0, 1]),
                title="methylation" if i == 4 else None,
                axis=alt.Axis(
                    labels=show_x,
                    ticks=True,
                    domain=True,
                ),
            ),
            y=alt.Y(
                "type:N",
                title="cell fate" if show_y else None,
                axis=alt.Axis(
                    labels=show_y,
                    ticks=False,
                    domain=True,
                ),
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
            width=200,
            height=150,
            title={
                "text": biomarker,
                "subtitle": score_map.get(biomarker, ""),
                "subtitleFontSize": 9,
            },
        )
    )

    charts.append(c)

chart = alt.concat(*charts, columns=ncols)
chart.save(snakemake.output[0])

print(f"Chart saved to {snakemake.output[0]}")
