import pandas as pd
import altair as alt


pd.set_option("display.max_columns", None)

df = pd.read_parquet(snakemake.input, engine="pyarrow")


# cg21699252 - 2:15938891-15938892 (-1)
# cg00933813 - 10:35594676-35594677 (-1)
# cg00661673 - 4:168841140-168841141 (1)
psc_position_pairs = [("2", 15938891), ("10", 35594676), ("4", 168841140)]

filtered_df = df[
    df.apply(
        lambda row: (row["chromosome"], row["position"]) in psc_position_pairs,
        axis=1,
    )
].reset_index(drop=True)

print(filtered_df)
# Initialisiere das Dictionary für die neuen Daten
pluripotency_scores = []

# Berechnungen für jede Gruppe (Biomarker)


# Erstellen der Zeilen für den aktuellen Biomarker
for methylation_type, methylation_col in [
    ("PSC", "psc_methylation"),
    ("MESO", "mesoderm_methylation"),
    ("ENDO", "endoderm_methylation"),
    ("ECTO", "ectoderm_methylation"),
    # ("ENDOMESO", "endomeso_methylation"),
]:
    pluripotency_scores.append(
        {
            "pluripotency_score": filtered_df[methylation_col].sum() / 100,
            "type": (
                "undifferentiated" if methylation_type == "PSC" else "differentiated"
            ),
            "line": methylation_type,
        }
    )

# Erstellen des neuen DataFrames
result_df = pd.DataFrame(pluripotency_scores)

plot_df = pd.DataFrame(pluripotency_scores)

# Altair-Plot erstellen
chart = (
    alt.Chart(plot_df)
    .mark_point(size=100)
    .encode(
        x=alt.X("type:N", title="Selection set"),
        y=alt.Y("pluripotency_score:Q", title="Score"),
        color=alt.Color(
            "line:N",
            title="Line",
            scale=alt.Scale(
                domain=["PSC", "MESO", "ENDO", "ECTO"],
                range=["black", "blue", "red", "green"],
            ),
        ),
        tooltip=["pluripotency_score", "type", "line"],  # Interaktive Tooltip-Infos
    )
    .properties(title="Pluripotency Scores by Type and Line", width=600, height=400)
)

chart.save(snakemake.output[0], scale_factor=2.0)
