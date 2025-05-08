rule df_from_calls:
    input:
        undifferentiated="results/{platform}/meth_calling/psc/calls.vcf",
        meso="results/{platform}/meth_calling/mesoderm/calls.vcf",
        endo="results/{platform}/meth_calling/endoderm/calls.vcf",
        ecto="results/{platform}/meth_calling/ectoderm/calls.vcf",
        # "results/{platform}/meth_calling/{group}/calls.vcf",
    output:
        "results/{platform}/meth_calling/calls.parquet",
    conda:
        "../envs/plot.yaml"
    script:
        "../scripts/df_from_calls.py"


rule most_variable_positions:
    input:
        "results/{platform}/meth_calling/calls.parquet",
    output:
        "results/{platform}/meth_calling/most_variable_positions.parquet",
    conda:
        "../envs/plot.yaml"
    script:
        "../scripts/most_variable.py"


rule scatter_plot:
    input:
        calls="results/{platform}/meth_calling/calls.parquet",
        # most_variable_positions="results/{platform}/meth_calling/most_variable_positions.parquet",
    output:
        report(
            "results/{platform}/plots_paper/{group2}/scatter_plot.html",
            caption="../report/scatter_plot.rst",
            category="Plots paper",
            subcategory="Plot 1B",
            labels=lambda wildcards: {
                "Platform": wildcards.platform,
                "Type": wildcards.group2,
            },
        ),
    params:
        group1=config["ref_sample"],
        group2=lambda wildcards: wildcards.group2,
    wildcard_constraints:
        group2="(?!endo_meso).*",
    conda:
        "../envs/plot.yaml"
    script:
        "../scripts/scatter_plot.py"


rule scatter_plot_pb_np:
    input:
        pacbio="results/nanopore/meth_calling/calls.parquet",
        nanopore="results/pacbio/meth_calling/calls.parquet",
    output:
        "results/comp_pb_np/meth_comp_pb_np_{group}.html",
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
        calls="results/{platform}/meth_calling/calls.parquet",
        most_variable_positions="results/{platform}/meth_calling/most_variable_positions.parquet",
    output:
        report(
            "results/{platform}/plots_paper/endo_meso/scatter_plot.html",
            caption="../report/scatter_plot.rst",
            category="Plots paper",
            subcategory="Plot 1C",
            labels=lambda wildcards: {
                "Platform": wildcards.platform,
                "Type": "endo_meso",
            },
        ),
    params:
        group1="mesoderm",
        group2="endoderm",
    conda:
        "../envs/plot.yaml"
    script:
        "../scripts/scatter_plot.py"


rule pluripotency_score_psc:
    input:
        "results/{platform}/meth_calling/calls.parquet",
    output:
        report(
            "results/{platform}/plots_paper/pluripotency_score_psc.html",
            caption="../report/scatter_plot.rst",
            category="Plots paper",
            subcategory="Plot 2B",
            labels=lambda wildcards: {
                "Platform": wildcards.platform,
                "Type": "undifferentiated",
            },
        ),
    conda:
        "../envs/plot.yaml"
    script:
        "../scripts/pluripotency_score_psc.py"


rule pluripotency_score_all:
    input:
        "results/{platform}/meth_calling/calls.parquet",
    output:
        report(
            "results/{platform}/plots_paper/pluripotency_score_all.html",
            caption="../report/scatter_plot.rst",
            category="Plots paper",
            subcategory="Plot 3B",
            labels=lambda wildcards: {
                "Platform": wildcards.platform,
                "Type": "differentiated",
            },
        ),
    conda:
        "../envs/plot.yaml"
    script:
        "../scripts/pluripotency_score_all_heatmap.py"
