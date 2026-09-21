rule plot_expression_to_meth:
    input:
        expr="results/tables/tpm-matrix/just_get_counts.tpm-matrix.sorted.tsv",
        meth="results/{platform}/{caller}/meth_calling/chipseeker_{fdr}.tsv",
    output:
        report(
            "results/{platform}/{caller}/plots_paper/expression_to_meth_{annotation}_{fdr}.pdf",
            caption="../report/expression_to_meth.rst",
            category="Expression vs. methylation",
            subcategory=lambda wildcards: wildcards.platform,
            labels=lambda wildcards: {
                "region": wildcards.annotation,
                "FDR": wildcards.fdr,
            },
        ),
    conda:
        "../envs/python.yaml"
    log:
        "logs/plot_expression_to_meth/{platform}_{caller}_{annotation}_{fdr}.log",
    script:
        "../scripts/plot_expression_to_meth.py"
