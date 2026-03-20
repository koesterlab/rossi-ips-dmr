import sys

import altair as alt
import pandas as pd

# Setup logging
sys.stderr = open(snakemake.log[0], "w", buffering=1)

# Configure pandas display
pd.set_option("display.max_columns", None)

# Read input data
df = pd.read_parquet(snakemake.input, engine="pyarrow")

# Define methylation column mapping
METHYLATION_COLS = {
    "psc_methylation": "psc",
    "endoderm_methylation": "endoderm",
    "mesoderm_methylation": "mesoderm",
    "ectoderm_methylation": "ectoderm",
}

# Define biomarker positions as DataFrame
BIOMARKER_POSITIONS = pd.DataFrame(
    [
        # PSC biomarkers
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
        # Endoderm biomarkers
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
        # Mesoderm biomarkers
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
        # Ectoderm biomarkers
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
        # Endomesoderm biomarkers
        {
            "chromosome": "5",
            "position": 111309143,
            "direction": 1,
            "cg_id": "cg23385847",
            "biomarker": "endomesoderm",
        },
        {
            "chromosome": "5",
            "position": 24208809,
            "direction": 1,
            "cg_id": "cg24919344",
            "biomarker": "endomesoderm",
        },
        {
            "chromosome": "10",
            "position": 33773306,
            "direction": 1,
            "cg_id": "cg11147278",
            "biomarker": "endomesoderm",
        },
    ]
)

# Ensure correct data types
BIOMARKER_POSITIONS["chromosome"] = BIOMARKER_POSITIONS["chromosome"].astype(str)

# Prepare data: add position_pair for filtering
df["chromosome"] = df["chromosome"].astype(str)
df["position_pair"] = list(zip(df["chromosome"], df["position"]))

biomarker_pairs = set(
    zip(BIOMARKER_POSITIONS["chromosome"], BIOMARKER_POSITIONS["position"])
)
filtered_df = df[df["position_pair"].isin(biomarker_pairs)].copy()

# Merge biomarker metadata
filtered_df = filtered_df.merge(
    BIOMARKER_POSITIONS[["chromosome", "position", "direction", "biomarker"]],
    on=["chromosome", "position"],
    how="left",
)

# Reshape to long format
long_df = filtered_df.melt(
    id_vars=["chromosome", "position", "direction", "biomarker"],
    value_vars=list(METHYLATION_COLS.keys()),
    var_name="layer",
    value_name="methylation",
).assign(layer=lambda x: x["layer"].map(METHYLATION_COLS))

print("Long format data:")
print(long_df)
print("\n")


# Helper function to calculate adjusted methylation score
def calculate_adjusted_methylation(methylation_value, direction):
    """
    Adjust methylation based on direction:
    - If direction=1: return methylation as is
    - If direction=-1: return 1 - methylation
    """
    return methylation_value if direction == 1 else 1 - methylation_value


# Define biomarker focus sets
BIOMARKER_FOCUS_SETS = {
    "psc": {"psc"},
    "endoderm": {"endoderm"},
    "mesoderm": {"mesoderm"},
    "ectoderm": {"ectoderm"},
    "endomesoderm": {"endoderm", "mesoderm"},
}

# Generate charts for each biomarker
charts = []
for biomarker_name in BIOMARKER_POSITIONS["biomarker"].unique():
    # Filter data for this biomarker
    biomarker_df = long_df[long_df["biomarker"] == biomarker_name].copy()

    # Determine focus type based on biomarker
    focus_set = BIOMARKER_FOCUS_SETS[biomarker_name]
    biomarker_df["type"] = biomarker_df["layer"].apply(
        lambda layer: "focus" if layer in focus_set else "unfocused"
    )

    # Calculate adjusted methylation scores
    biomarker_df["adjusted_methylation"] = biomarker_df.apply(
        lambda row: calculate_adjusted_methylation(
            row["methylation"], row["direction"]
        ),
        axis=1,
    )

    # Print statistics
    print(f"Biomarker: {biomarker_name}")
    print(biomarker_df)

    score_per_layer = biomarker_df.groupby("layer")["adjusted_methylation"].sum()
    print(f"Score per layer (sum of adjusted methylation):")
    print(score_per_layer)
    print("\n")

    # Create visualization
    chart = (
        alt.Chart(biomarker_df)
        .mark_point(size=100)
        .encode(
            x=alt.X("methylation:Q", title="Methylation"),
            y=alt.Y("type:N", title="Selection set"),
            color=alt.Color("layer:N", title="Germ Layer"),
        )
        .properties(
            title=f"Methylation Scores by Layer - {biomarker_name}",
            width=100,
            height=150,
        )
    )
    charts.append(chart)

# Save concatenated charts
alt.hconcat(*charts).save(snakemake.output[0])
print(f"Chart saved to {snakemake.output[0]}")
