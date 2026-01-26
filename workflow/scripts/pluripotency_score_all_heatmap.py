import pandas as pd
import altair as alt

sys.stderr = open(snakemake.log[0], "w", buffering=1)


pd.set_option("display.max_columns", None)


df = pd.read_parquet(snakemake.input, engine="pyarrow")
methylation_cols = {
    "psc_methylation": "psc",
    "endoderm_methylation": "endoderm",
    "mesoderm_methylation": "mesoderm",
    "ectoderm_methylation": "ectoderm",
}

psc_position_pairs = {
    ("2", 15938891): "cg21699252",
    ("10", 35594676): "cg00933813",
    ("4", 168841140): "cg00661673",
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
# for positions in [
#     ("PSC", psc_position_pairs),
#     ("ENDO", endo_position_pairs),
#     ("MESO", meso_position_pairs),
#     ("ECTO", ecto_position_pairs),
#     ("ENDOMESO", endomeso_position_pairs),
# ]:
#     filtered_df = df[
#         df.apply(
#             lambda row: (row["chromosome"], row["position"]) in psc_position_pairs,
#             axis=1,
#         )
#     ].reset_index(drop=True)
#     meth_vals = []

#     long_df = filtered_df.melt(
#         id_vars=["chromosome", "position"],
#         value_vars=list(methylation_cols.keys()),
#         var_name="type",
#         value_name="methylation",
#     ).assign(type=lambda x: x["type"].map(methylation_cols))
#     print(long_df)
#     # df["position_pair"] = list(zip(df["chromosome"].astype(str), df["position"]))
#     # filtered_df = df[df["position_pair"].isin(positions[1])].reset_index(drop=True)
#     # filtered_df["cg_id"] = filtered_df["position_pair"].map(positions[1])
#     # print(positions[0])
#     # print(filtered_df)
#     # print("\n\n")

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

long_df = filtered_df.melt(
    id_vars=["chromosome", "position"],
    value_vars=list(methylation_cols.keys()),
    var_name="layer",
    value_name="methylation",
).assign(layer=lambda x: x["layer"].map(methylation_cols))
print(long_df)


charts = []
for biomarker, position_pairs in [
    ("psc", psc_position_pairs),
    ("endoderm", endo_position_pairs),
    ("mesoderm", meso_position_pairs),
    ("ectoderm", ecto_position_pairs),
    ("endomesoderm", endomeso_position_pairs),
]:
    biomarker_df = long_df[
        long_df.apply(
            lambda row: (row["chromosome"], row["position"]) in position_pairs,
            axis=1,
        )
    ].reset_index(drop=True)
    biomarker_df["type"] = biomarker_df.apply(
        lambda row: (
            "focus"
            if (
                row["layer"] == biomarker
                or (
                    biomarker == "endomesoderm"
                    and row["layer"] in {"endoderm", "mesoderm"}
                )
            )
            else "unfocused"
        ),
        axis=1,
    )
    print(f"Biomarker: {biomarker}")
    print(biomarker_df)
    print("\n\n")

    chart = (
        alt.Chart(biomarker_df)
        .mark_point(size=100)
        .encode(
            x=alt.X("type:N", title="Selection set"),
            y=alt.Y("methylation:Q", title="Methylation"),
            color=alt.Color(
                "layer:N",
                title="Germ Layer",
                # scale=alt.Scale(
                #     domain=["PSC", "MESO", "ENDO", "ECTO"],
                #     range=["black", "blue", "red", "green"],
                # ),
            ),
        )
        .properties(
            title=f"Methylation Scores by Layer {biomarker}", width=600, height=400
        )
    )
    charts.append(chart)
alt.vconcat(*charts).save(snakemake.output[0], scale_factor=2.0)
# df_melted = filtered_df.melt(
#     id_vars=["position", "chromosome", "cg_id"],
#     value_vars=[
#         "psc_methylation",
#         "endoderm_methylation",
#         "mesoderm_methylation",
#         "ectoderm_methylation",
#     ],
#     var_name="Cell_Type",
#     value_name="Methylation",
# )
# df_melted["Cell_Type"] = df_melted["Cell_Type"].replace(
#     {
#         "psc_methylation": "PSC",
#         "endoderm_methylation": "Endoderm",
#         "mesoderm_methylation": "Mesoderm",
#         "ectoderm_methylation": "Ectoderm",
#     }
# )
# df_melted["biomarker"] = df_melted.apply(
#     lambda row: (
#         "PSC"
#         if (row["chromosome"], row["position"]) in psc_position_pairs.keys()
#         else (
#             "ENDO"
#             if (row["chromosome"], row["position"]) in endo_position_pairs.keys()
#             else (
#                 "MESO"
#                 if (row["chromosome"], row["position"]) in meso_position_pairs.keys()
#                 else (
#                     "ECTO"
#                     if (row["chromosome"], row["position"])
#                     in ecto_position_pairs.keys()
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
# psc_values = df_melted[df_melted["Cell_Type"] == "PSC"][
#     ["position", "Methylation"]
# ].set_index("position")["Methylation"]
# df_melted["Methylation_diff"] = df_melted.apply(
#     lambda row: row["Methylation"] - psc_values.get(row["position"], 0), axis=1
# )
# charts = []
# for biomarker in ["PSC", "ENDO", "MESO", "ENDOMESO", "ECTO"]:
#     df = df_melted[df_melted["biomarker"] == biomarker]
#     heatmap = (
#         alt.Chart(df)
#         .mark_rect()
#         .encode(
#             x=alt.X("Cell_Type:N", title="Cell Type"),
#             y=alt.Y(
#                 "cg_id:N",
#                 title="Genomic Position",
#             ),
#             color=alt.Color(
#                 "Methylation_diff:Q", scale=alt.Scale(scheme="redblue", domain=[-1, 1])
#             ),
#             tooltip=[
#                 "position",
#                 "chromosome",
#                 "Cell_Type",
#                 "Methylation",
#                 "Methylation_diff",
#                 "biomarker",
#             ],
#         )
#         .properties(width=500, height=400, title=f"Methylation Heatmap {biomarker}")
#     )
#     charts.append(heatmap)
# alt.vconcat(*charts).save(snakemake.output[0], scale_factor=2.0)
