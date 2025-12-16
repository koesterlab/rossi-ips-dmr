from collections import defaultdict
import pandas as pd
import sys

sys.stderr = open(snakemake.log[0], "w", buffering=1)

pd.set_option("display.max_columns", None)

df = pd.read_parquet(snakemake.input[0], engine="pyarrow")
group2 = snakemake.params["group2"]

with open(str(snakemake.output[0]), "w") as outfile:
    outfile.write(f"chr\tpos\tpsc\t{group2}\n")
    for _, row in df.iterrows():
        chrom = row["chromosome"]
        pos = row["position"]
        methylation_psc = row["psc_methylation"] / 100.0
        methylation = row[f"{group2}_methylation"] / 100.0
        outfile.write(f"{chrom}\t{pos}\t")
        outfile.write(f"{methylation_psc}\t")
        outfile.write(f"{methylation}\t")
        outfile.write("\n")
