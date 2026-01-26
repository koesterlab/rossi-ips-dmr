import pandas as pd
import re
import sys
import pysam

sys.stderr = open(snakemake.log[0], "w", buffering=1)


def compute_bias(format_values):
    """Returns the bias label if present in FORMAT values, else 'normal'."""
    bias_labels = ["SB", "ROB", "RPB", "SCB", "HE", "ALB"]
    if any(value != "." for value in format_values[6:12]):
        return bias_labels[format_values[6:12].index(".")]
    return "normal"


def read_tool_file(file_path, axis_name, meth_caller="varlo"):
    df = []

    if meth_caller == "varlo":
        vcf = pysam.VariantFile(file_path)

        for record in vcf:
            chrom = record.chrom
            position = record.pos
            alt = record.alts[0] if record.alts else None

            if alt != "<METH>":
                continue

            sample_afs = []
            sample_dps = []
            sample_bias = "normal"

            for sample in record.samples.values():
                try:
                    af = float(sample.get("AF", [0])[0])
                    dp = int(sample.get("DP", 0))
                    bias = compute_bias([v[0] if isinstance(v, tuple) and len(v) == 1 else v for v in sample.values()])
                    if bias != "normal":
                        sample_bias = bias
                    sample_afs.append(af * 100)
                    sample_dps.append(dp)
                except Exception:
                    pass

            meth_rate = sum(sample_afs) / len(sample_afs) if sample_afs else 0
            coverage = sum(sample_dps) / len(sample_dps) if sample_dps else 0

            info = record.info
            alpha = float(snakemake.params["alpha"])

            def phred_to_prob(score):
                return 10 ** (-float(score) / 10)

            if "PROB_PRESENT" in info:
                if isinstance(info["PROB_PRESENT"][0], float):
                    prob_present = phred_to_prob(info["PROB_PRESENT"][0])
                else:
                    # print(f"Missing probability info at {chrom}:{position}", file=sys.stderr)   
                    continue
            else:
                # print(f"Missing probability info at {chrom}:{position}", file=sys.stderr)
                continue

            prob_absent = phred_to_prob(info["PROB_ABSENT"][0])
            prob_artifact = phred_to_prob(info["PROB_ARTIFACT"][0])
            prob_identified = max(prob_present, prob_absent + prob_artifact)

            if prob_identified < (1 - alpha):
                continue

            df.append(
                [chrom, position, meth_rate, coverage, sample_bias, prob_identified]
            )

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
        f"{axis_name}_prob_identified",
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

print("Starting to read files")
undiff_df = read_tool_file(undiff_file, "psc", meth_caller)
print("Finished going through undiff")
meso_df = read_tool_file(meso_file, "mesoderm", meth_caller)
print("Finished going through meso")
endo_df = read_tool_file(endo_file, "endoderm", meth_caller)
print("Finished going through endo")
ecto_df = read_tool_file(ecto_file, "ectoderm", meth_caller)
print("Finished going through ecto")


df = merge_dfs(undiff_df, meso_df, endo_df, ecto_df)
print("Finished merging", file=sys.stderr)


df.to_parquet(snakemake.output[0], engine="pyarrow", compression="snappy")
