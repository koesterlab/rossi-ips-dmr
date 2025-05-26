

rule download_metilene:
    output:
        temp("resources/tools/metilene_v02-8.tar.gz"),
    shell:
        """
        wget -O {output} http://www.bioinf.uni-leipzig.de/Software/metilene/metilene_v02-8.tar.gz
        """


rule unpack_metilene:
    input:
        "resources/tools/metilene_v02-8.tar.gz",
    output:
        directory("resources/tools/metilene"),
    shell:
        """
        mkdir -p {output}
        tar -xzf {input} -C {output} --strip-components=1
        """


rule metilene_input:
    input:
        "results/{platform}/{caller}/meth_calling/calls.parquet",
    output:
        "results/{platform}/{caller}/dmr_calls/{group2}/metilene_input.txt",
    params:
        group2=lambda wildcards: wildcards.group2,
    conda:
        "../envs/plot.yaml"
    log:
        "logs/metilene_input/{platform}/{caller}/{group2}.log",
    script:
        "../scripts/metilene_input.py"


rule call_metilene:
    input:
        "results/{platform}/{caller}/dmr_calls/{group2}/metilene_input.txt",
    output:
        "results/{platform}/{caller}/dmr_calls/{group2}/metilene_output.bed",
    conda:
        "../envs/metilene.yaml"
    threads: 10
    params:
        base_exp_number=config["ref_sample"],
    shell:
        """
        metilene -t {threads} -c 2 -a {params.base_exp_number} -b {wildcards.group2} {input} > {output}
        """


rule plot_dmrs:
    input:
        met=directory("resources/tools/metilene"),
        met_out="results/{platform}/{caller}/dmr_calls/{group2}/metilene_output.bed",
    output:
        "results/{platform}/{caller}/dmr_calls/{group2}/plots/dmr_qval.0.05.bedgraph",
        report(
            "results/{platform}/{caller}/dmr_calls/{group2}/plots/dmr_qval.0.05.pdf",
            caption="../report/metilene_plots.rst",
            category="DMR plots",
            subcategory=lambda wildcards: f"Metilene: {wildcards.platform} - {wildcards.caller}",
            labels=lambda wildcards: {
                "experiment 1": config["ref_sample"],
                "experiment 2": wildcards.group2,
            },
        ),
    conda:
        "../envs/metilene.yaml"
    params:
        path_prefix=lambda wildcards: f"results/{wildcards.platform}/{wildcards.caller}/dmr_calls/{wildcards.group2}/plots/dmr",
        sample=lambda wildcards: {wildcards.group2.split("-")[0]},
        base_exp_number=config["ref_sample"],
    shell:
        """
        PARENT_DIR=$(dirname {output})/dmr
        echo {params.base_exp_number}
        echo $PARENT_DIR
        perl {input.met}/metilene_output.pl -q {input.met_out} -o $PARENT_DIR -a {params.base_exp_number} -b {wildcards.group2}
        """


# Still included for debugging reason but not used right now
rule _plot_pvals:
    input:
        met_out="results/{platform}/{caller}/dmr_calls/{group2}/metilene_output.bed",
    output:
        "results/{platform}/{caller}/dmr_calls/{group2}/plots/pvals.html",
    conda:
        "../envs/plot.yaml"
    params:
        path_prefix=lambda wildcards: f"results/{platform}/{caller}/dmr_calls/{wildcards.group2}",
        sample=lambda wildcards: {wildcards.group2.split("-")[0]},
    script:
        "../scripts/plot_pvals.py"
