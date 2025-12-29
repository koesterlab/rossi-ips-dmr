import polars as pl
import numpy as np

sample_rows = []
unit_rows = []

for fastq in snakemake.input.fastqs:
    sample_name = fastq.split("/")[-1].split(".")[0]
    cell_type = snakemake.params.types[sample_name]
    sample_rows.append({"sample": f"{cell_type}_{sample_name[-2:]}", "condition": cell_type, "batch_effect": "NA"})

    unit_rows.append(
        {
            "sample": f"{cell_type}_{sample_name[-2:]}",
            "unit": 1,
            "fragment_len_mean": "NA",
            "fragment_len_sd": "NA",
            "fq1": fastq,
            "fq2": "NA",
            "bam_single": "NA",
            "bam_paired": "NA",
            "fastp_adapters": "NA",
            "fastp_extra": "NA",
        }
    )


samples_df = pl.DataFrame(sample_rows)
units_df = pl.DataFrame(unit_rows)

samples_df.write_csv(snakemake.output["samples"], separator="\t")
units_df.write_csv(snakemake.output["units"], separator="\t")
