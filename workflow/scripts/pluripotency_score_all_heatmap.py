import pandas as pd
import altair as alt


pd.set_option("display.max_columns", None)


df = pd.read_parquet(snakemake.input, engine="pyarrow")


# cg21699252 - 2:15938891-15938892 (-1)
# cg00933813 - 10:35594676-35594677 (-1)
# cg00661673 - 4:168841140-168841141 (1)

# Die chromosomen-Positionen, nach denen gefiltert werden soll
psc_position_pairs = [("2", 15938891), ("10", 35594676), ("4", 168841140)]
endo_position_pairs = [("6", 12886978), ("11", 8840472), ("8", 125637563)]
meso_position_pairs = [("2", 128638846), ("17", 15966293), ("12", 122872581)]
endomeso_position_pairs = [("5", 111309143), ("5", 24208809), ("10", 33773306)]
ecto_position_pairs = [("15", 71314056), ("5", 107661548), ("14", 68569167)]


all_positions = (
    psc_position_pairs
    + endo_position_pairs
    + meso_position_pairs
    + ecto_position_pairs
    + endomeso_position_pairs
)
print(all_positions)
# Filter anwenden
filtered_df = df[
    df.apply(
        lambda row: (row["chromosome"], row["position"]) in all_positions,
        axis=1,
    )
].reset_index(drop=True)

# filtered_df["biomarker"] = filtered_df.apply(
#     lambda row: (
#         "PSC"
#         if (row["chromosome"], row["position"]) in psc_position_pairs
#         else (
#             "ENDO"
#             if (row["chromosome"], row["position"]) in endo_position_pairs
#             else (
#                 "MESO"
#                 if (row["chromosome"], row["position"]) in meso_position_pairs
#                 else (
#                     "ECTO"
#                     if (row["chromosome"], row["position"]) in ecto_position_pairs
#                     else (
#                         "ENDOMESO"
#                         if (row["chromosome"], row["position"])
#                         in endomeso_position_pairs
#                         else None
#                     )
#                 )
#             )
#         )
#     ),
#     axis=1,
# )

# print(filtered_df)

df_melted = filtered_df.melt(
    id_vars=["position", "chromosome"],
    value_vars=[
        "psc_methylation",
        "endoderm_methylation",
        "mesoderm_methylation",
        "ectoderm_methylation",
    ],
    var_name="Cell_Type",
    value_name="Methylation",
)

# Zelltypen schöner benennen
df_melted["Cell_Type"] = df_melted["Cell_Type"].replace(
    {
        "psc_methylation": "PSC",
        "endoderm_methylation": "Endoderm",
        "mesoderm_methylation": "Mesoderm",
        "ectoderm_methylation": "Ectoderm",
    }
)

print(df_melted)

df_melted["biomarker"] = df_melted.apply(
    lambda row: (
        "PSC"
        if (row["chromosome"], row["position"]) in psc_position_pairs
        else (
            "ENDO"
            if (row["chromosome"], row["position"]) in endo_position_pairs
            else (
                "MESO"
                if (row["chromosome"], row["position"]) in meso_position_pairs
                else (
                    "ECTO"
                    if (row["chromosome"], row["position"]) in ecto_position_pairs
                    else (
                        "ENDOMESO"
                        if (row["chromosome"], row["position"])
                        in endomeso_position_pairs
                        else None
                    )
                )
            )
        )
    ),
    axis=1,
)

psc_values = df_melted[df_melted["Cell_Type"] == "PSC"][
    ["position", "Methylation"]
].set_index("position")["Methylation"]

# Nun für jede Zeile die Differenz berechnen, indem der Methylierungswert von der PSC Position subtrahiert wird
df_melted["Methylation_diff"] = df_melted.apply(
    lambda row: row["Methylation"] - psc_values.get(row["position"], 0), axis=1
)


print(df_melted)

# Altair Heatmap
charts = []
for biomarker in ["PSC", "ENDO", "MESO", "ENDOMESO", "ECTO"]:
    df = df_melted[df_melted["biomarker"] == biomarker]
    heatmap = (
        alt.Chart(df)
        .mark_rect()
        .encode(
            x=alt.X("Cell_Type:N", title="Cell Type"),
            y=alt.Y(
                "position:N",
                title="Genomic Position",
                # sort=alt.EncodingSortField(field="position", op="count", order="descending"),
            ),
            color=alt.Color("Methylation_diff:Q", scale=alt.Scale(scheme="redblue")),
            tooltip=[
                "position",
                "Cell_Type",
                "Methylation",
                "Methylation_diff",
                "biomarker",
            ],
        )
        .properties(width=500, height=400, title=f"Methylation Heatmap {biomarker}")
    )
    charts.append(heatmap)
alt.vconcat(*charts).save(snakemake.output[0], scale_factor=2.0)
# chart = (
#     alt.Chart(plot_df)
#     .mark_point(size=100)
#     .encode(
#         x=alt.X("biomarker:N", title="Biomarker"),
#         y=alt.Y("pluripotency_score:Q", title="Score"),
#         color=alt.Color(
#             "line:N",
#             title="Line",
#             scale=alt.Scale(
#                 domain=["PSC", "MESO", "ENDO", "ECTO", "ENDOMESO"],
#                 range=["black", "blue", "red", "green", "orange"],
#             ),
#         ),
#         tooltip=["pluripotency_score", "type", "line"],
#     )
#     .properties(title="Pluripotency Scores by Type and Line", width=600, height=400)
# )
# biomarkers = ["PSC", "MESO", "ENDO", "ECTO", "ENDOMESO"]
# pluripotency_scores = []
# for biomarker in biomarkers:
#     biomarker_df = filtered_df[filtered_df["biomarker"] == biomarker]
#     print("biomarker_df", biomarker_df)
#     for methylation_type, methylation_col in [
#         ("PSC", "psc_methylation"),
#         ("MESO", "mesoderm_methylation"),
#         ("ENDO", "endoderm_methylation"),
#         ("ECTO", "ectoderm_methylation"),
#         ("ENDOMESO", "endomeso_methylation"),
#     ]:
#         pluripotency_scores.append(
#             {
#                 "pluripotency_score": biomarker_df[methylation_col].sum() / 100,
#                 "type": (
#                     "undifferentiated"
#                     if methylation_type == "PSC"
#                     else "differentiated"
#                 ),
#                 "line": methylation_type,
#                 "biomarker": biomarker,
#             }
#         )
# plot_df = pd.DataFrame(pluripotency_scores)
# print(plot_df)
# chart = (
#     alt.Chart(plot_df)
#     .mark_point(size=100)
#     .encode(
#         x=alt.X("biomarker:N", title="Biomarker"),
#         y=alt.Y("pluripotency_score:Q", title="Score"),
#         color=alt.Color(
#             "line:N",
#             title="Line",
#             scale=alt.Scale(
#                 domain=["PSC", "MESO", "ENDO", "ECTO", "ENDOMESO"],
#                 range=["black", "blue", "red", "green", "orange"],
#             ),
#         ),
#         tooltip=["pluripotency_score", "type", "line"],
#     )
#     .properties(title="Pluripotency Scores by Type and Line", width=600, height=400)
# )
# chart.save(snakemake.output[0], scale_factor=2.0)
