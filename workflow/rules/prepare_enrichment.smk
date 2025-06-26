rule get_ensembl_genes:
    input:
        "results/{platform}/{caller}/dmr_calls/{group2}/genes_transcripts/chipseeker_complete.tsv",
    output:
        "results/{platform}/{caller}/dmr_calls/{group2}/genes_transcripts/ensembl_genes.tsv",
    conda:
        "../envs/biomart.yaml"
    params:
        species=get_bioc_species_name(),
        version=config["resources"]["ref"]["release"],
    log:
        "logs/get_ensembl_genes/{platform}/{caller}/dmr_calls/{group2}/genes_transcripts.log",
    script:
        "../scripts/get_ensembl_genes.R"


# rule compose_sample_sheet:
#     input:
#         config["samples"],
#         config["units"],
#         kallisto_output=kallisto_output,
#     output:
#         "results/{platform}/{caller}/sleuth/{group}.samples.tsv",
#     log:
#         "logs/{group}.compose-sample-sheet.log",
#     params:
#         units=units,
#         samples=samples,
#     group:
#         "sleuth-init"
#     script:
#         "../scripts/compose-sample-sheet.py"
