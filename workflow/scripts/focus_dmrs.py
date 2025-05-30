import pybedtools
import sys

sys.stderr = open(snakemake.log[0], "w")

ecto = pybedtools.BedTool(snakemake.input["ecto"])
endo = pybedtools.BedTool(snakemake.input["endo"])
meso = pybedtools.BedTool(snakemake.input["meso"])

ecto_only = ecto.intersect(endo, v=True).intersect(meso, v=True)
endo_only = endo.intersect(ecto, v=True).intersect(meso, v=True)
meso_only = meso.intersect(ecto, v=True).intersect(endo, v=True)

exclusive_sets = [
    (ecto_only, snakemake.output["ecto"]),
    (endo_only, snakemake.output["endo"]),
    (meso_only, snakemake.output["meso"]),
]

for regions, outfile in exclusive_sets:
    with open(outfile, "w") as out:
        for region in regions:
            try:
                value = float(region[4])

                if abs(value) > 0.20:
                    print(region, file=out, end="")
            except (IndexError, ValueError):
                continue
