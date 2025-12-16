# Define (un)methylated and prob_present as in compute_precision_recall

import pandas as pd
import sys


# Redirect standard error to snakemake log file
sys.stderr = open(snakemake.log[0], "w", buffering=1)


# Set up plotting and display
pd.set_option("display.max_columns", None)

# Threshold for methylation (When do we say a site is methylated. 0 is bad since the truth (avg bedgraph) has nearly never meth rates of 0)
methylation_threshold = snakemake.params["meth_threshold"]

df = pd.read_parquet(snakemake.input, engine="pyarrow")

lineages = ["psc", "ectoderm", "mesoderm", "endoderm"]

for method in lineages:
    # Todo: Think about this again. How to look at prob_presents
    df[f"{method}_meth_flag"] = df[f"{method}_methylation"].apply(
        lambda x: 1 if x > methylation_threshold else 0
    )


# Keep rows where exactly one meth_flag is 1 or exactly  one is 0
meth_flag_cols = [f"{method}_meth_flag" for method in lineages]
df = df[(df[meth_flag_cols].sum(axis=1) == 1) | (df[meth_flag_cols].sum(axis=1) == 3)]
# print(df, file=sys.stderr)
print(df, file=sys.stderr)
df.to_parquet(snakemake.output[0], engine="pyarrow", compression="snappy")
