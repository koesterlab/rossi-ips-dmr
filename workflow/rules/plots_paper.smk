rule pluripotency_score_all:
    input:
        "results/{platform}/{caller}/meth_calling/calls_1.0.parquet",
    output:
        report(
            "results/{platform}/{caller}/plots_paper/pluripotency_score_all.{plot_type}",
            caption="../report/scatter_plot.rst",
            category="Plots paper",
            subcategory=lambda wildcards: f"{wildcards.platform} - {wildcards.caller}",
            labels=lambda wildcards: {
                "Plot": "3B",
                "Type": "differentiated",
            },
        ),
    conda:
        "../envs/python.yaml"
    resources:
        mem_mb=16000,
    log:
        "logs/pluripotency_score_all/{platform}_{caller}_{plot_type}.log",
    script:
        "../scripts/pluripotency_score_all_heatmap.py"


# Compares DMRs between PacBio and Nanopore
rule dmr_scatter_comparison:
    input:
        pacbio=expand(
            "results/pacbio/varlo/base_psc/dmr_calls/{group2}/genes_transcripts/1.0/chipseeker_postprocessed.tsv",
            group2=get_non_base_layers("psc"),
        ),
        nanopore=expand(
            "results/nanopore/varlo/base_psc/dmr_calls/{group2}/genes_transcripts/1.0/chipseeker_postprocessed.tsv",
            group2=get_non_base_layers("psc"),
        ),
    output:
        report(
            "results/platforms_combined/varlo/plots_paper/scatter_comparison.{plot_type}",
            caption="../report/scatter_comparison.rst",
            category="PacBio vs. Nanopore",
            labels=lambda wildcards: {
                "Type": "platform comparison",
            },
        ),
    conda:
        "../envs/python.yaml"
    log:
        "logs/dmr_scatter_comparison/{plot_type}.log",
    resources:
        mem_mb=16000,
    script:
        "../scripts/dmr-heatmap_comparison.py"


rule concatenate_figure:
    input:
        plots=[
            "results/platforms_combined/varlo/plots_paper/scatter_comparison.pdf",
            "results/platforms_combined/varlo/base_psc/rna_new/diffexp_vs_dmrs_promoter.pdf",
            "results/platforms_combined/varlo/base_psc/rna_new/diffexp_vs_dmrs_unfiltered.pdf",
        ],
    output:
        "results/platforms_combined/varlo/plots_paper/concatenated.{plot_type}",
    conda:
        "../envs/fitz.yaml"
    log:
        "logs/concatenate_figure/{plot_type}.log",
    resources:
        mem_mb=4000,
    script:
        "../scripts/concat_plots.py"
