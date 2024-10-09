rule get_ensembl_genes:
    input:
        "results/dmr_calls/{group2}/genes_transcripts/chipseeker.tsv",
    output:
        "results/dmr_calls/{group2}/genes_transcripts/ensembl_genes.tsv",
    conda:
        "../envs/biomart.yaml"
    log:
        "logs/get_ensembl_genes_{group2}.log",
    params:
        species=get_bioc_species_name(),
        version=config["resources"]["ref"]["release"],
    script:
        "../scripts/get_ensembl_genes.R"


# We only merge with polars because the join in get_ensembl_genes.R does not work
rule annotate_chipseeker:
    input:
        chipseeker="results/dmr_calls/{group2}/genes_transcripts/chipseeker.tsv",
        genes="results/dmr_calls/{group2}/genes_transcripts/ensembl_genes.tsv",
    output:
        "results/dmr_calls/{group2}/genes_transcripts/chipseeker_postprocessed.tsv",
    conda:
        "../envs/python_standard.yaml"
    log:
        "logs/annotate_chipseeker{group2}.log",
    script:
        "../scripts/annotate_chipseeker.py"


# rule compose_sample_sheet:
#     input:
#         config["samples"],
#         config["units"],
#         kallisto_output=kallisto_output,
#     output:
#         "results/sleuth/{group}.samples.tsv",
#     log:
#         "logs/{group}.compose-sample-sheet.log",
#     params:
#         units=units,
#         samples=samples,
#     group:
#         "sleuth-init"
#     script:
#         "../scripts/compose-sample-sheet.py"
