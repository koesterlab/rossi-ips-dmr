import sys

import polars as pl

sys.stderr = open(snakemake.log[0], "w", buffering=1)
pl.Config.set_tbl_rows(10000)
pl.Config.set_tbl_cols(300)
nanopore = pl.read_parquet(snakemake.input.nanopore)
print(nanopore, file="/mnt/workspace/rossi-ips-dmr/test_nanopore.csv")
pacbio = pl.read_parquet(snakemake.input.pacbio)
print(pacbio, file="/mnt/workspace/rossi-ips-dmr/test_pacbio.csv")

methylation_cols = [col for col in nanopore.columns if col.endswith("_methylation")]
print(nanopore, file=sys.stderr)


def filter_df(df: pl.DataFrame, methylation_cols: list[str]) -> pl.DataFrame:
    """Keep only rows where exactly one sample is fully methylated (>95)
    and all others are fully unmethylated (<5), or vice versa."""
    n = len(methylation_cols)
    return (
        df.with_columns(
            count_high=pl.sum_horizontal([pl.col(c) > 95 for c in methylation_cols]),
            count_low=pl.sum_horizontal([pl.col(c) < 5 for c in methylation_cols]),
        )
        .filter(
            ((pl.col("count_high") == 1) & (pl.col("count_low") == n - 1))
            | ((pl.col("count_high") == n - 1) & (pl.col("count_low") == 1))
        )
        .drop(["count_high", "count_low"])
    )


nanopore_filtered = filter_df(nanopore, methylation_cols)
pacbio_filtered = filter_df(pacbio, methylation_cols)

# Inner join on chromosome + position. Explicit suffix so we know exactly
# which columns come from pacbio_filtered, regardless of which columns
# happen to collide (methylation_cols always will; coverage_cols only if
# nanopore and pacbio coverage frames share column names).
df_combined = nanopore_filtered.join(
    pacbio_filtered, on=["chromosome", "position"], how="inner", suffix="_right"
)

# Keep only positions where the methylation call agrees between nanopore and pacbio.
# Agreement = both above 95 ("methylated") or both below 5 ("unmethylated").
for col in methylation_cols:
    right_col = f"{col}_right"
    df_combined = df_combined.filter(
        ((pl.col(col) > 55) & (pl.col(right_col) > 60))
        | ((pl.col(col) < 45) & (pl.col(right_col) < 45))
    ).drop(right_col)

df_combined = df_combined.select(
    pl.col("chromosome"),
    pl.col("position"),
    *methylation_cols,
)
print(df_combined, file=sys.stderr)
df_combined.write_csv(snakemake.output[0])
