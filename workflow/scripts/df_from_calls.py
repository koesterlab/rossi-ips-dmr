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

                chrom, position, alternative, info_field, format_field, values = (
                    parts[0],
                    int(parts[1]),
                    parts[4],
                    parts[7],
                    parts[8],
                    parts[9].split(":"),
                )

                if alternative == "<METH>":
                    format_fields = format_field.split(":")
                    dp_index = format_fields.index("DP")
                    af_index = format_fields.index("AF")

                    meth_rate = float(values[af_index]) * 100
                    coverage = int(values[dp_index])

                    match = re.search(r"PROB_PRESENT=([\d\.]+)", info_field)
                    prob_present = 0
                    try:
                        prob_present = 10 ** (-float(match.group(1)) / 10)
                    except:
                        print("Prob present not found")

                    bias = compute_bias(values)

                    df.append([chrom, position, meth_rate, coverage, bias, prob_present])
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
        f"{axis_name}_prob_present",
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
