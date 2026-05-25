ruleorder:
    scatter_plot_endo_meso > scatter_plot

rule scatter_plot:
    input:
        calls="results/{platform}/{caller}/meth_calling/calls.parquet",
    output:
        report(
            "results/{platform}/{caller}/base_{base}/plots_paper/{group2}/scatter_plot.png",
            caption="../report/scatter_plot.rst",
            category="Plots paper",
            subcategory=lambda wildcards: f"{wildcards.platform} - {wildcards.caller}",
            labels=lambda wildcards: {
                "Plot": "1B",
                "Base": wildcards.base,
                "Type": wildcards.group2,
            },
        ),
    params:
        group1=lambda wildcards: wildcards.base,
        group2=lambda wildcards: wildcards.group2,
        meth_caller=lambda wildcards: wildcards.caller,
    wildcard_constraints:
        group2="(?!endo_meso).*",
    resources:
        mem_mb=16000,
    conda:
        "../envs/python.yaml"
    log:
        "logs/scatter_plot/{platform}_{caller}_{base}_{group2}.log",
    script:
        "../scripts/scatter_plot.py"


rule scatter_plot_endo_meso:
    input:
        calls="results/{platform}/{caller}/meth_calling/calls.parquet",
    output:
        report(
            "results/{platform}/{caller}/plots_paper/endo_meso/scatter_plot.png",
            caption="../report/scatter_plot.rst",
            category="Plots paper",
            subcategory=lambda wildcards: f"{wildcards.platform} - {wildcards.caller}",
            labels=lambda wildcards: {
                "Plot": "1C",
                "Type": "endo_meso",
            },
        ),
    resources:
        mem_mb=16000,
    params:
        group1="mesoderm",
        group2="endoderm",
        meth_caller=lambda wildcards: wildcards.caller,
    conda:
        "../envs/python.yaml"
    log:
        "logs/scatter_plot_endo_meso/{platform}_{caller}.log",
    script:
        "../scripts/scatter_plot.py"


rule pluripotency_score_all:
    input:
        "results/{platform}/{caller}/meth_calling/calls.parquet",
    output:
        report(
            "results/{platform}/{caller}/plots_paper/pluripotency_score_all.html",
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
        "logs/pluripotency_score_all/{platform}_{caller}.log",
    script:
        "../scripts/pluripotency_score_all_heatmap.py"


rule dmr_heatmap_comparison:
    input:
        pacbio = lambda wildcards: expand(
            "results/pacbio/varlo/base_psc/dmr_calls/{group2}/genes_transcripts/chipseeker_postprocessed.tsv",
            group2=get_non_base_layers("psc"),
        ),
        nanopore = lambda wildcards: expand(
            "results/nanopore/varlo/base_psc/dmr_calls/{group2}/genes_transcripts/chipseeker_postprocessed.tsv",
            group2=get_non_base_layers("psc"),
        ),
    output:
        # report(
        "results/platforms_combined/varlo/plots_paper/heatmaps_comparison.png",
        #     caption="../report/heatmap.rst",
        #     category="DMR plots",
        #     subcategory=lambda wildcards: f"Heatmaps: {wildcards.platform} - {wildcards.caller}",
        #     labels=lambda wildcards: {
        #         "base": wildcards.base,
        #         "genetic element": wildcards.type,
        #     },
        # ),
    conda:
        "../envs/plot.yaml"
    log:
        "logs/dmr_heatmap/comparison.log",
    resources:
        mem_mb=16000,
    params:
        base = "psc"
    script:
        "../scripts/dmr-heatmap_comparison.py"
