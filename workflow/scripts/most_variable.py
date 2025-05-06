import pandas as pd
import numpy as np


df = pd.read_parquet(snakemake.input, engine="pyarrow")


methylation_cols = [
    "psc_methylation",
    "mesoderm_methylation",
    "endoderm_methylation",
    "ectoderm_methylation",
]

df_methylation = df[["chromosome", "position"] + methylation_cols]


df_methylation["std_dev"] = df_methylation[methylation_cols].std(axis=1)


top_10000 = df_methylation.nlargest(10000, "std_dev").reset_index(drop=True)

top_10000.to_parquet(snakemake.output[0], engine="pyarrow", compression="snappy")
