rule plot_expression_to_meth:
    input:
        expr="results/tables/tpm-matrix/just_get_counts_psc.tpm-matrix.sorted.tsv",
        meth="results/{platform}/{caller}/meth_calling/chipseeker_{fdr}.tsv"
    output:
        "results/{platform}/{caller}/plots_paper/expression_to_meth_{annotation}_{fdr}.pdf"
    conda:
        "../envs/python.yaml"
    log:
        "logs/plot_expression_to_meth/{platform}_{caller}_{annotation}_{fdr}.log"
    params:
        annotation=lambda wildcards: wildcards.annotation
    script:
        "../scripts/plot_expression_to_meth.py"
