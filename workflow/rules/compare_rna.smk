rule filter_rna_seq:
    input:
        rna_seq="resources/rna_seq/rna_seq.xlsx",
    output:
        rna_seq_filtered="resources/rna_seq/rna_seq_filtered_computed.xlsx",
    conda:
        "../envs/filter_rna.yaml"
    log:
        "logs/filter_rna_seq.log",
    script:
        "../scripts/filter_rna_seq.py"


rule rna_seq:
    input:
        genes_transcripts="results/{platform}/{caller}/dmr_calls/{group2}/genes_transcripts/chipseeker_postprocessed_complete.tsv",
        # rna_seq="resources/rna_seq/rna_seq_filtered.xlsx",
        rna_seq="resources/rna_seq/rna_seq_filtered_computed.xlsx",
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
        barplot=report(
            "results/{platform}/{caller}/rna_seq/{group2}_barplot.html",
            caption="../report/annotations.rst",
            htmlindex="index.html",
            category="RNA Seq Comparison",
            subcategory=lambda wildcards: f"Barplot: {wildcards.platform} - {wildcards.caller}",
            labels=lambda wildcards: {
                "experiment 1": config["ref_sample"],
                "experiment 2": wildcards.group2,
            },
        ),
    params:
        germ_layer=lambda wildcards: wildcards.group2,
    conda:
        "../envs/plot_rna.yaml"
    log:
        "logs/rna_seq/{platform}_{caller}_{group2}.log",
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
    log:
        "logs/rna_seq_annotations/{platform}_{caller}_{group2}.log",
    wrapper:
        "v3.13.0/utils/datavzrd"
