import pybedtools
import sys

sys.stderr = open(snakemake.log[0], "w")

ecto = pybedtools.BedTool(snakemake.input["ecto"])
endo = pybedtools.BedTool(snakemake.input["endo"])
meso = pybedtools.BedTool(snakemake.input["meso"])
all_files = [ecto, endo, meso]

common_regions = ecto.intersect(endo).intersect(meso)

# Remember regions to avoid duplicates
seen = [set(), set(), set()]

with open(snakemake.output["ecto"], "w") as out_ecto, open(
    snakemake.output["endo"], "w"
) as out_endo, open(snakemake.output["meso"], "w") as out_meso:

    outputs = [out_ecto, out_endo, out_meso]

    for region in common_regions:
        region_bt = pybedtools.BedTool([region])

        overlaps = [f.intersect(region_bt, wa=True) for f in all_files]

        try:
            meth_values = [float(overlap[0].fields[4]) for overlap in overlaps]
        except (IndexError, ValueError):
            continue

        diff_meth = sum(abs(val) > 0.05 for val in meth_values)

        if diff_meth in (1, 3):
            for i, (match, seen_set, out) in enumerate(zip(overlaps, seen, outputs)):
                line = str(match[0])
                if line not in seen_set:
                    print(line, file=out, end="")
                    seen_set.add(line)
