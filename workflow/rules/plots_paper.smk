rule df_from_calls:
    input:
        undifferentiated="results/{platform}/{caller}/meth_calling/psc/{caller}.vcf",
        meso="results/{platform}/{caller}/meth_calling/mesoderm/{caller}.vcf",
        endo="results/{platform}/{caller}/meth_calling/endoderm/{caller}.vcf",
        ecto="results/{platform}/{caller}/meth_calling/ectoderm/{caller}.vcf",
        # "results/{platform}/{caller}/meth_calling/{group}/calls.vcf",
    output:
        "results/{platform}/{caller}/meth_calling/calls.parquet",
    conda:
        "../envs/plot.yaml"
    log:
        "logs/df_from_calls/{platform}/{caller}.log"
    params:
        meth_caller=lambda wildcards: wildcards.caller,
        prob_pres_threshhold=config["prob_pres_threshold"],
        prob_abs_threshhold=config["prob_abs_threshold"],
    script:
        "../scripts/df_from_calls.py"

rule focus_df_on_diff_methylated_loci:
    input:
        "results/{platform}/{caller}/meth_calling/calls.parquet",
    output:
        "results/{platform}/{caller}/meth_calling/calls_focused.parquet",
    conda:
        "../envs/plot.yaml"
    log:
        "logs/focus_df_on_diff_methylated_loci/{platform}/{caller}.log"
    params:
        meth_threshold=config["meth_threshold"],
    script:
        "../scripts/focus_loci.py"


rule most_variable_positions:
    input:
        "results/{platform}/{caller}/meth_calling/calls_focused.parquet",
    output:
        "results/{platform}/{caller}/meth_calling/most_variable_positions.parquet",
    conda:
        "../envs/plot.yaml"
    script:
        "../scripts/most_variable.py"


rule scatter_plot:
    input:
        calls="results/{platform}/{caller}/meth_calling/calls_focused.parquet",
        # most_variable_positions="results/{platform}/{caller}/meth_calling/most_variable_positions.parquet",
    output:
        report(
            "results/{platform}/{caller}/plots_paper/{group2}/scatter_plot.html",
            caption="../report/scatter_plot.rst",
            category="Plots paper",
            subcategory=lambda wildcards: f"{wildcards.platform} - {wildcards.caller}",
            labels=lambda wildcards: {
                "Plot": "1B",
                "Type": wildcards.group2,
            },
        ),
    params:
        group1=config["ref_sample"],
        group2=lambda wildcards: wildcards.group2,
        meth_caller=lambda wildcards: wildcards.caller
    wildcard_constraints:
        group2="(?!endo_meso).*",
    conda:
        "../envs/plot.yaml"
    script:
        "../scripts/scatter_plot.py"


rule scatter_plot_pb_np:
    input:
        pacbio="results/nanopore/{caller}/meth_calling/calls_focused.parquet",
        nanopore="results/pacbio/{caller}/meth_calling/calls_focused.parquet",
    output:
        "results/comp_pb_np/meth_comp_pb_np_{group}.png",
    params:
        group=lambda wildcards: wildcards.group,
    wildcard_constraints:
        # group2="(?!endo_meso).*",
    conda:
        "../envs/plot.yaml"
    script:
        "../scripts/scatter_plot_pb_np.py"


rule scatter_plot_endo_meso:
    input:
        # Must be list because uses same plot function as scatter_plot
        calls="results/{platform}/{caller}/meth_calling/calls_focused.parquet",
        # most_variable_positions="results/{platform}/{caller}/meth_calling/most_variable_positions.parquet",
    output:
        report(
            "results/{platform}/{caller}/plots_paper/endo_meso/scatter_plot.html",
            caption="../report/scatter_plot.rst",
            category="Plots paper",
            subcategory=lambda wildcards: f"{wildcards.platform} - {wildcards.caller}",
            labels=lambda wildcards: {
                "Plot": "1C",
                "Type": "endo_meso",
            },
        ),
    params:
        group1="mesoderm",
        group2="endoderm",
        meth_caller=lambda wildcards: wildcards.caller
    conda:
        "../envs/plot.yaml"
    script:
        "../scripts/scatter_plot.py"


rule pluripotency_score_psc:
    input:
        "results/{platform}/{caller}/meth_calling/calls_focused.parquet",
    output:
        report(
            "results/{platform}/{caller}/plots_paper/pluripotency_score_psc.html",
            caption="../report/scatter_plot.rst",
            category="Plots paper",
            subcategory=lambda wildcards: f"{wildcards.platform} - {wildcards.caller}",
            labels=lambda wildcards: {
                "Plot": "2B",
                "Type": "undifferentiated",
            },
        ),
    conda:
        "../envs/plot.yaml"
    script:
        "../scripts/pluripotency_score_psc.py"


rule pluripotency_score_all:
    input:
        "results/{platform}/{caller}/meth_calling/calls_focused.parquet",
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
        "../envs/plot.yaml"
    script:
        "../scripts/pluripotency_score_all_heatmap.py"
