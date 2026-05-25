import sys

import polars as pl

sys.stderr = snakemake.log[0]

coverage = pl.empty()
for cov_file in snakemake.input:
    df = (
        pl.read_csv(
            cov_file,
            separator="\t",
            has_header=False,
            new_columns=["chrom", "start", "end", "cov"],
        )
        .with_columns((pl.col("start") + pl.col("end")) / 2)
        .alias("pos")
    )

    coverage = (
        pl.merge(coverage, df, on=["chrom", "pos"], how="outer")
        .fill_null(0)
        .select(["chrom", "pos", "cov"])
    )

coverage.to_parquet(snakemake.output[0])
