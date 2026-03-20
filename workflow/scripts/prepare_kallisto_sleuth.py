import sys

import polars as pl

sys.stderr = open(snakemake.log[0], "w", buffering=1)

sample_rows = []
unit_rows = []
fastqs_old = snakemake.input.old_fastqs
fastqs_new = snakemake.input.new_fastqs
fastqs = fastqs_old + fastqs_new

# Define contrast columns and their mappings
contrasts = {
    "ectoderm_vs_psc_old": ("ectoderm", "psc", "old"),
    "endoderm_vs_psc_old": ("endoderm", "psc", "old"),
    "mesoderm_vs_psc_old": ("mesoderm", "psc", "old"),
    "ectoderm_vs_endoderm_old": ("ectoderm", "endoderm", "old"),
    "ectoderm_vs_mesoderm_old": ("ectoderm", "mesoderm", "old"),
    "mesoderm_vs_endoderm_old": ("mesoderm", "endoderm", "old"),
    "ectoderm_vs_psc_new": ("ectoderm", "psc", "new"),
    "endoderm_vs_psc_new": ("endoderm", "psc", "new"),
    "mesoderm_vs_psc_new": ("mesoderm", "psc", "new"),
    "ectoderm_vs_endoderm_new": ("ectoderm", "endoderm", "new"),
    "ectoderm_vs_mesoderm_new": ("ectoderm", "mesoderm", "new"),
    "mesoderm_vs_endoderm_new": ("mesoderm", "endoderm", "new"),
}

for fastq in fastqs:
    sample_name = fastq.split("/")[-1].split(".")[0]
    relative_path = fastq.split("resources/")[1]
    cell_type = (
        snakemake.params.labels_old[sample_name]
        if fastq in fastqs_old
        else snakemake.params.labels_new[sample_name]
    )
    rna_data = "old" if fastq in fastqs_old else "new"

    # Build contrast columns
    contrast_values = {}
    for col_name, (pos_type, neg_type, data_type) in contrasts.items():
        if rna_data == data_type:
            if cell_type == pos_type:
                contrast_values[col_name] = "+"
            elif cell_type == neg_type:
                contrast_values[col_name] = "-"
            else:
                contrast_values[col_name] = ""
        else:
            contrast_values[col_name] = ""

    sample_rows.append(
        {
            "sample": f"{cell_type}_{sample_name[-2:]}",
            "label": cell_type,
            "rna_data": rna_data,
            **contrast_values,
        }
    )

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
