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
            x=alt.X("methylation:Q", title="Methylation"),
            y=alt.Y("type:N", title="Selection set"),
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
            title=f"Methylation Scores by Layer {biomarker}", width=100, height=150
        )
    )
    charts.append(chart)
alt.hconcat(*charts).save(snakemake.output[0])
