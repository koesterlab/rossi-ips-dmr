import pandas as pd
import re


def compute_bias(prob, format_values):
    probs = re.split("[=;]", prob)
    probs_values = [probs[i] for i in [1, 3, 5]]
    try:
        [float(value) for value in probs_values]
    except:
        return "bias"

    # min_value = min(prob_values)
    # min_index = prob_values.index(min_value)

    bias_labels = ["SB", "ROB", "RPB", "SCB", "HE", "ALB"]

    if any(value != "." for value in format_values[6:12]):
        bias = bias_labels[format_values[6:12].index(".")]
    else:
        bias = "normal"

    # if min_index != 0:
    #     bias = "normal"
    return bias


def read_call_file(file_path, axis_name):
    df = []
    with open(file_path, "r") as tool_file:

        for line in tool_file:
            if line.startswith("#") or line.startswith("track"):
                continue
            parts = line.strip().split("\t")

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
                prob_present = 1000
                try:
                    prob_present = 10 ** (-float(match.group(1)) / 10)
                except:
                    print("Prob present not found")

                bias = compute_bias(info_field, values)

                df.append([chrom, position, meth_rate, coverage, bias, prob_present])
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

undiff_file = snakemake.input["undifferentiated"]
meso_file = snakemake.input["meso"]
endo_file = snakemake.input["endo"]
ecto_file = snakemake.input["ecto"]

undiff_df = read_call_file(undiff_file, "psc")
meso_df = read_call_file(meso_file, "mesoderm")
endo_df = read_call_file(endo_file, "endoderm")
ecto_df = read_call_file(ecto_file, "ectoderm")

print(undiff_df)
print(meso_df)
print(endo_df)
print(ecto_df)


df = merge_dfs(undiff_df, meso_df, endo_df, ecto_df)

print(df)
df.to_parquet(snakemake.output[0], engine="pyarrow", compression="snappy")
