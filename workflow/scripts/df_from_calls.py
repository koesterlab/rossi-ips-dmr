import pandas as pd
import re
import sys

sys.stderr = open(snakemake.log[0], "w")


def compute_bias(format_values):
    """Returns the bias label if present in FORMAT values, else 'normal'."""
    bias_labels = ["SB", "ROB", "RPB", "SCB", "HE", "ALB"]
    if any(value != "." for value in format_values[6:12]):
        return bias_labels[format_values[6:12].index(".")]
    return "normal"


def read_tool_file(file_path, axis_name, meth_caller="varlo"):
    df = []
    with open(file_path, "r") as tool_file:

        for line in tool_file:
            if line.startswith("#") or line.startswith("track"):
                continue
            parts = line.strip().split("\t")

            if meth_caller == "varlo":



                chrom = parts[0]
                position = int(parts[1])
                alt = parts[4]
                info_field = parts[7]
                format_field = parts[8]

                if alt != "<METH>":
                    continue

                # Parse sample AF values
                format_fields = format_field.split(":")
                try:
                    af_index = format_fields.index("AF")
                    dp_index = format_fields.index("DP")
                except ValueError:
                    continue

                sample_afs = []
                sample_dps = []
                sample_bias = "normal"
                for sample_fmt in parts[9:]:
                    values = sample_fmt.split(":")
                    try:
                        af = float(values[af_index])
                        dp = int(values[dp_index])
                        bias = compute_bias(values)
                        if bias != "normal":
                            sample_bias = bias

                        sample_afs.append(af * 100)
                        sample_dps.append(dp)
                    except (ValueError, IndexError):
                        pass

                meth_rate = sum(sample_afs) / len(sample_afs) if sample_afs else 0
                coverage = sum(sample_dps) / len(sample_dps) if sample_dps else 0

                # Parse INFO fields
                info_dict = dict(
                    item.split("=", 1)
                    for item in info_field.strip().split(";")
                    if "=" in item
                )

                alpha = float(snakemake.params["alpha"])

                def phred_to_prob(score):
                    return 10 ** (-float(score) / 10)

                def is_missing(val):
                    return val is None or val == "."

                # Determine probability of methylation
                if is_missing(info_dict.get("PROB_PRESENT")):
                    if is_missing(info_dict.get("PROB_HIGH")) or is_missing(
                        info_dict.get("PROB_LOW")
                    ):
                        print(
                            f"Missing probability info at {chrom}:{position}",
                            file=sys.stderr,
                        )
                        continue
                    prob_present = phred_to_prob(
                        info_dict["PROB_HIGH"]
                    ) + phred_to_prob(info_dict["PROB_LOW"])
                else:
                    prob_present = phred_to_prob(info_dict["PROB_PRESENT"])

                prob_absent = phred_to_prob(info_dict.get("PROB_ABSENT", 0))
                prob_artifact = phred_to_prob(info_dict.get("PROB_ARTIFACT", 0))

                if max(prob_present, prob_absent + prob_artifact) < (1 - alpha):
                    print(
                        f"Low confidence site skipped: {chrom}:{position}",
                        file=sys.stderr,
                    )
                    continue

                df.append([chrom, position, meth_rate, coverage, sample_bias])
            elif meth_caller == "modkit":
                chrom = parts[0].removeprefix("chr")
                position = int(parts[2])
                details = parts[9].split()
                coverage = int(details[0])
                meth_rate = float(details[1])
                df.append(
                    [
                        chrom,
                        position,
                        meth_rate,
                        coverage,
                        "normal",
                        0,
                    ]
                )
            elif meth_caller == "pb_CpG_tools":
                chrom = parts[0]
                position = int(parts[2])
                meth_rate = float(parts[3])
                coverage = int(parts[5])
                df.append(
                    [
                        chrom,
                        position,
                        meth_rate,
                        coverage,
                        "normal",
                        0,
                    ]
                )
    columns = [
        "chromosome",
        "position",
        f"{axis_name}_methylation",
        f"{axis_name}_coverage",
        f"{axis_name}_bias",
        # f"{axis_name}_prob_present",
    ]
    return pd.DataFrame(df, columns=columns)


def merge_dfs(undifferentiated_df, meso_df, endo_df, ecto_df):
    undifferentiated_df_filtered = undifferentiated_df[
        undifferentiated_df["psc_bias"] == "normal"
    ]
    meso_df_filtered = meso_df[meso_df["mesoderm_bias"] == "normal"]
    endo_df_filtered = endo_df[endo_df["endoderm_bias"] == "normal"]
    ecto_df_filtered = ecto_df[ecto_df["ectoderm_bias"] == "normal"]

    merged_df = (
        undifferentiated_df_filtered.merge(
            endo_df_filtered, on=["chromosome", "position"], how="outer"
        )
        .merge(meso_df_filtered, on=["chromosome", "position"], how="outer")
        .merge(ecto_df_filtered, on=["chromosome", "position"], how="outer")
    )

    return merged_df


pd.set_option("display.max_columns", None)
meth_caller = snakemake.params["meth_caller"]
undiff_file = snakemake.input["undifferentiated"]
meso_file = snakemake.input["meso"]
endo_file = snakemake.input["endo"]
ecto_file = snakemake.input["ecto"]

undiff_df = read_tool_file(undiff_file, "psc", meth_caller)
print("Finished going through undiff", file=sys.stderr)
meso_df = read_tool_file(meso_file, "mesoderm", meth_caller)
print("Finished going through meso", file=sys.stderr)
endo_df = read_tool_file(endo_file, "endoderm", meth_caller)
print("Finished going through endo", file=sys.stderr)
ecto_df = read_tool_file(ecto_file, "ectoderm", meth_caller)
print("Finished going through ecto", file=sys.stderr)


df = merge_dfs(undiff_df, meso_df, endo_df, ecto_df)
print("Finished merging", file=sys.stderr)


df.to_parquet(snakemake.output[0], engine="pyarrow", compression="snappy")
