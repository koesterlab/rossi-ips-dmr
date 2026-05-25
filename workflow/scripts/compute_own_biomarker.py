import sys

import polars as pl

sys.stderr = open(snakemake.log[0], "w", buffering=1)

pl.Config.set_tbl_rows(300)
pl.Config.set_tbl_cols(300)
nanopore = pl.read_parquet(snakemake.input.nanopore)
pacbio = pl.read_parquet(snakemake.input.pacbio)
coverage_nanopore = pl.read_parquet(snakemake.input.coverage_nanopore)
coverage_pacbio = pl.read_parquet(snakemake.input.coverage_pacbio)

nanopore = pl.merge(
    nanopore, coverage_nanopore, on=["chromosome", "position"], how="left"
)
pacbio = pl.merge(pacbio, coverage_pacbio, on=["chromosome", "position"], how="left")

methylation_cols = [col for col in nanopore.columns if col.endswith("_methylation")]
coverage_cols = [col for col in nanopore.columns if col.endswith("_coverage")]
coverage_cols_right = [
    col + "_right" for col in nanopore.columns if col.endswith("_coverage")
]


# Filter beide
def filter_df(df, methylation_cols):
    return (
        df.with_columns(
            count_100=pl.concat_list(methylation_cols).list.count_matches(100),
            count_0=pl.concat_list(methylation_cols).list.count_matches(0),
        )
        .filter(
            (
                (pl.col("count_100") == 1)
                & (pl.col("count_0") == len(methylation_cols) - 1)
            )
            | (
                (pl.col("count_100") == len(methylation_cols) - 1)
                & (pl.col("count_0") == 1)
            )
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

df_combined = df_combined.select(
    pl.col("chromosome"),
    pl.col("position"),
    *methylation_cols,
    *coverage_cols,
    *coverage_cols_right,
)
print(df_combined, file=sys.stderr)
df_combined.write_parquet(snakemake.output[0])
