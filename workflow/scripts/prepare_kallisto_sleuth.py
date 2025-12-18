import polars as pl
import numpy as np

sample_rows = []
unit_rows = []

for fastq in snakemake.input.fastqs:
    sample_name = fastq.split("/")[-1].split(".")[0]
    cell_type = snakemake.params.types[sample_name]
    sample_rows.append({"sample": sample_name, "condition": cell_type})

    unit_rows.append(
        {
            "sample": sample_name,
            "unit": 1,
            "fragment_len_mean": "?",
            "fragment_len_sd": "?",
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
