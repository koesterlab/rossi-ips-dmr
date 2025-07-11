rule rna_seq:
    input:
        genes_transcripts="results/{platform}/{caller}/dmr_calls/{group2}/genes_transcripts/chipseeker_postprocessed_complete.tsv",
        rna_seq="resources/rna_seq/rna_seq.xlsx",
    output:
        table="results/{platform}/{caller}/rna_seq/{group2}.tsv",
        scatter=report(
            "results/{platform}/{caller}/rna_seq/{group2}_scatter.html",
            caption="../report/annotations.rst",
            htmlindex="index.html",
            category="RNA Seq Comparison",
            subcategory=lambda wildcards: f"Scatter: {wildcards.platform} - {wildcards.caller}",
            labels=lambda wildcards: {
                "experiment 1": config["ref_sample"],
                "experiment 2": wildcards.group2,
            },
        ),
        histogram=report(
            "results/{platform}/{caller}/rna_seq/{group2}_histogram.html",
            caption="../report/annotations.rst",
            htmlindex="index.html",
            category="RNA Seq Comparison",
            subcategory=lambda wildcards: f"Scatter: {wildcards.platform} - {wildcards.caller}",
            labels=lambda wildcards: {
                "experiment 1": config["ref_sample"],
                "experiment 2": wildcards.group2,
            },
        ),
    params:
        germ_layer=lambda wildcards: wildcards.group2,
    conda:
        "../envs/plot.yaml"
    log:
        "logs/rna_seq/{platform}/{caller}/{group2}.log",
    script:
        "../scripts/compare_rna_seq.py"


rule rna_seq_annotations:
    input:
        config=workflow.source_path("../resources/rna_annotated.yaml"),
        rna_annotations="results/{platform}/{caller}/rna_seq/{group2}.tsv",
    output:
        report(
            directory("results/{platform}/{caller}/rna_seq/datavzrd/{group2}"),
            caption="../report/annotations.rst",
            htmlindex="index.html",
            category="RNA Seq Comparison",
            subcategory=lambda wildcards: f"Tables: {wildcards.platform} - {wildcards.caller}",
            labels=lambda wildcards: {
                "experiment 1": config["ref_sample"],
                "experiment 2": wildcards.group2,
                "type": "table",
            },
        ),
    params:
        base_experiment=lambda wildcards: config["ref_sample"],
    wrapper:
        "v3.13.0/utils/datavzrd"
