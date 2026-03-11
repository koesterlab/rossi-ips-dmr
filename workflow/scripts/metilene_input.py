import sys

import pandas as pd

sys.stderr = open(snakemake.log[0], "w", buffering=1)

pd.set_option("display.max_columns", None)

df = pd.read_parquet(snakemake.input[0], engine="pyarrow")
base = snakemake.params["base"]
group2 = snakemake.params["group2"]

with open(str(snakemake.output[0]), "w") as outfile:
    outfile.write(f"chr\tpos\t{base}\t{group2}\n")
    for _, row in df.iterrows():
        chrom = row["chromosome"]
        pos = row["position"]
        methylation_base = row[f"{base}_methylation"] / 100.0
        methylation_group2 = row[f"{group2}_methylation"] / 100.0
        outfile.write(f"{chrom}\t{pos}\t")
        outfile.write(f"{methylation_base}\t")
        outfile.write(f"{methylation_group2}\t")
        outfile.write("\n")
