import sys

import polars as pl

sys.stderr = open(snakemake.log[0], "w", buffering=1)
pl.Config.set_tbl_rows(100)
pl.Config.set_tbl_cols(-1)

# A site counts as fully (un)methylated in a germ layer above/below these values (%)
FULLY_METHYLATED = 95
FULLY_UNMETHYLATED = 5
# Thresholds for the agreement between Nanopore and PacBio (%)
AGREE_METHYLATED = 55
AGREE_UNMETHYLATED = 45

nanopore = pl.read_parquet(snakemake.input.nanopore)
pacbio = pl.read_parquet(snakemake.input.pacbio)
methylation_cols = [col for col in nanopore.columns if col.endswith("_methylation")]


def filter_specific_sites(df: pl.DataFrame) -> pl.DataFrame:
    """Keep only sites where exactly one germ layer is fully methylated and all
    others are fully unmethylated, or vice versa."""
    n = len(methylation_cols)
    count_high = pl.sum_horizontal([pl.col(c) > FULLY_METHYLATED for c in methylation_cols])
    count_low = pl.sum_horizontal([pl.col(c) < FULLY_UNMETHYLATED for c in methylation_cols])
    return df.filter(
        ((count_high == 1) & (count_low == n - 1))
        | ((count_high == n - 1) & (count_low == 1))
    )


df_combined = filter_specific_sites(nanopore).join(
    filter_specific_sites(pacbio),
    on=["chromosome", "position"],
    how="inner",
    suffix="_pacbio",
)

# Keep only sites where Nanopore and PacBio agree in every germ layer.
df_combined = df_combined.filter(
    pl.all_horizontal(
        ((pl.col(col) > AGREE_METHYLATED) & (pl.col(f"{col}_pacbio") > AGREE_METHYLATED))
        | ((pl.col(col) < AGREE_UNMETHYLATED) & (pl.col(f"{col}_pacbio") < AGREE_UNMETHYLATED))
        for col in methylation_cols
    )
).select("chromosome", "position", *methylation_cols)

print(df_combined, file=sys.stderr)
df_combined.write_csv(snakemake.output[0])
