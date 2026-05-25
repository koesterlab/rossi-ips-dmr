import polars as pl

nanopore = pl.read_parquet(snakemake.input.nanopore)
pacbio = pl.read_parquet(snakemake.input.pacbio)

methylation_cols = [col for col in nanopore.columns if col.endswith("_methylation")]


# Filter beide
def filter_df(df, methylation_cols):
    return (
        df.with_columns(
            count_100=pl.concat_list(methylation_cols).list.count_matches(100),
            count_0=pl.concat_list(methylation_cols).list.count_matches(0),
        )
        .filter(
            (pl.col("count_100") == 1)
            & (pl.col("count_0") == len(methylation_cols) - 1)
        )
        .drop(["count_100", "count_0"])
    )


nanopore_filtered = filter_df(nanopore, methylation_cols)
pacbio_filtered = filter_df(pacbio, methylation_cols)

# Inner join auf chromosome + position
df_combined = nanopore_filtered.join(
    pacbio_filtered, on=["chromosome", "position"], how="inner"
)

# Behalte nur Zeilen wo alle methylation_cols identisch sind
for col in methylation_cols:
    df_combined = df_combined.filter(pl.col(col) == pl.col(f"{col}_right"))
    df_combined = df_combined.drop(f"{col}_right")

print(df_combined)
