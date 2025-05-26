import pandas as pd
import altair as alt


pd.set_option("display.max_columns", None)


df = pd.read_parquet(snakemake.input, engine="pyarrow")

print(df)
# cg21699252 - 2:15938891-15938892 (-1)
# cg00933813 - 10:35594676-35594677 (-1)
# cg00661673 - 4:168841140-168841141 (1)

# Die chromosomen-Positionen, nach denen gefiltert werden soll
psc_position_pairs = {
    ("2", 15938891): "cg21699252",
    ("10", 35594676): "cg00933813",
    ("4", 168841140): "cg00661673",
    # usw. – ergänze die restlichen Mappings
}
endo_position_pairs = {
    ("6", 12886978): "cg20548013",
    ("11", 8840472): "cg14521421",
    ("8", 125637563): "cg08913523",
}
meso_position_pairs = {
    ("2", 128638846): "cg14708360",
    ("17", 15966293): "cg08826152",
    ("12", 122872581): "cg11599718",
}
endomeso_position_pairs = {
    ("5", 111309143): "cg23385847",
    ("5", 24208809): "cg24919344",
    ("10", 33773306): "cg11147278",
}
ecto_position_pairs = {
    ("15", 71314056): "cg01907071",
    ("5", 107661548): "cg18118164",
    ("14", 68569167): "cg13075942",
}

all_positions = (
    psc_position_pairs
    | endo_position_pairs
    | meso_position_pairs
    | ecto_position_pairs
    | endomeso_position_pairs
)


df["position_pair"] = list(zip(df["chromosome"].astype(str), df["position"]))
filtered_df = df[df["position_pair"].isin(all_positions)].reset_index(drop=True)
filtered_df["cg_id"] = filtered_df["position_pair"].map(all_positions)
print("filtered_df")
print(filtered_df)

df_melted = filtered_df.melt(
    id_vars=["position", "chromosome", "cg_id"],
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


df_melted["biomarker"] = df_melted.apply(
    lambda row: (
        "PSC"
        if (row["chromosome"], row["position"]) in psc_position_pairs.keys()
        else (
            "ENDO"
            if (row["chromosome"], row["position"]) in endo_position_pairs.keys()
            else (
                "MESO"
                if (row["chromosome"], row["position"]) in meso_position_pairs.keys()
                else (
                    "ECTO"
                    if (row["chromosome"], row["position"])
                    in ecto_position_pairs.keys()
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
                "cg_id:N",
                title="Genomic Position",
                # sort=alt.EncodingSortField(field="position", op="count", order="descending"),
            ),
            color=alt.Color(
                "Methylation_diff:Q", scale=alt.Scale(scheme="redblue", domain=[-1, 1])
            ),
            tooltip=[
                "position",
                "chromosome",
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
