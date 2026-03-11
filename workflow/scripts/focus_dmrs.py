import sys

import pybedtools

sys.stderr = open(snakemake.log[0], "w", buffering=1)

this = pybedtools.BedTool(snakemake.input["this"])
other1 = pybedtools.BedTool(snakemake.input["other1"])
other2 = pybedtools.BedTool(snakemake.input["other2"])

this_only = this.intersect(other1, v=True).intersect(other2, v=True)

with open(snakemake.output[0], "w") as out:
    for region in this_only:
        try:
            print(region, file=out, end="")
        except (IndexError, ValueError):
            continue
