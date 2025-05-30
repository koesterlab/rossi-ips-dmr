

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


# | chr | start | stop | q-value | mean methylation difference | #CpGs | p (MWU) | p (2D KS) | mean g1 | mean g2 |
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


rule focus_dmrs:
    input:
        ecto="results/{platform}/{caller}/dmr_calls/ectoderm/metilene_output.bed",
        meso="results/{platform}/{caller}/dmr_calls/endoderm/metilene_output.bed",
        endo="results/{platform}/{caller}/dmr_calls/mesoderm/metilene_output.bed",
    output:
        ecto="results/{platform}/{caller}/dmr_calls/ectoderm/metilene_output_focused.bed",
        meso="results/{platform}/{caller}/dmr_calls/endoderm/metilene_output_focused.bed",
        endo="results/{platform}/{caller}/dmr_calls/mesoderm/metilene_output_focused.bed",
    conda:
        "../envs/bedtools.yaml"
    log:
        "logs/focus_dmrs/{platform}/{caller}.log",
    params:
        meth_threshold=config["meth_threshold"],
    script:
        "../scripts/focus_dmrs.py"


rule plot_dmrs:
    input:
        met="resources/tools/metilene",
        met_out="results/{platform}/{caller}/dmr_calls/{group2}/metilene_output.bed",
    output:
        bed="results/{platform}/{caller}/dmr_calls/{group2}/plots/dmr_qval.0.05.bedgraph",
        pdf=report(
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
        PARENT_DIR=$(dirname {output.bed})/dmr
        echo {params.base_exp_number}
        echo $PARENT_DIR
        perl {input.met}/metilene_output.pl -q {input.met_out} -o $PARENT_DIR -a {params.base_exp_number} -b {wildcards.group2}
        """
