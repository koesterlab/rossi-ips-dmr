import polars as pl
import sys

# sys.stderr = open(snakemake.log[0], "w", buffering=1)

pl.Config.set_tbl_rows(50)
pl.Config.set_tbl_cols(100)

rna_seq_df = pl.read_excel(snakemake.input.rna_seq)

# pseudocount
pc = 0.001

rna_seq_df = rna_seq_df.rename(
    {
        "SYMBOL": "genes",
    }
)
rna_seq_df = rna_seq_df.drop_nulls(subset=["genes"])


for layer in ["ecto", "endo", "meso"]:
    rna_seq_df = rna_seq_df.with_columns(
        [
            ((pl.col(f"{layer}1") + pc) / (pl.col("pluri1") + pc))
            .log(2)
            .alias(f"log2fc_{layer}1"),
            ((pl.col(f"{layer}2") + pc) / (pl.col("pluri2") + pc))
            .log(2)
            .alias(f"log2fc_{layer}2"),
            (
                (
                    ((pl.col(f"{layer}1") + pl.col(f"{layer}2")) / 2 + pc)
                    / ((pl.col("pluri1") + pl.col("pluri2")) / 2 + pc)
                ).log(2)
                / 2
            ).alias(f"log2fc_{layer}_avg_direct"),
        ]
    )

    rna_seq_df = rna_seq_df.with_columns(
        [
            (((pl.col(f"log2fc_{layer}1") + pl.col(f"log2fc_{layer}2")) / 2) / 2).alias(
                f"log2fc_{layer}_avg_combined"
            ),
        ]
    )

    # print(
    #     rna_seq_df.sort(f"log2fc_{layer}_avg_combined", descending=True).select(
    #         "genes",
    #         "pluri1",
    #         "pluri2",
    #         f"{layer}1",
    #         f"{layer}2",
    #         f"log2fc_{layer}1",
    #         f"log2fc_{layer}2",
    #         f"log2fc_{layer}_avg_combined",
    #         f"log2fc_{layer}_avg_direct",
    #     )
    # )
print(rna_seq_df)
rna_seq_df.write_excel(snakemake.output[0])
