# TODO: Reactivate, right now it deletes too much data
rule filter_calls:
    input:
        "results/{sample}/calls.bcf",
    output:
        "results/{sample}/calls.filtered.bcf",
    log:
        "logs/filter_calls_{sample}.log",
    conda:
        "envs/varlociraptor.yaml"
    params:
        event="PRESENT",
    shell:
        """
        cd {params.varlo_path}
        cargo run --release -- filter-calls control-fdr --mode local-smart {input} --events {params.event} --fdr 0.01 > {output}
        """


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
        expand(
            "results/{base_experiment}/calls.vcf", base_experiment=config["ref_sample"]
        ),
        "results/{group2}/calls.vcf",
    output:
        "results/dmr_calls/{group2}/metilene_input.txt",
    log:
        "../logs/{group2}/metilene_input.log",
    script:
        "../scripts/metilene_input.py"


rule call_metilene:
    input:
        "results/dmr_calls/{group2}/metilene_input.txt",
    output:
        "results/dmr_calls/{group2}/metilene_output.bed",
    log:
        "../logs/{group2}/call_metilene.log",
    conda:
        "../envs/metilene.yaml"
    threads: 10
    params:
        base_exp_number=config["ref_sample_number"],
    shell:
        """
        metilene -t {threads} -c 2 -a {params.base_exp_number} -b {wildcards.group2} {input} > {output}
        """


rule plot_dmrs:
    input:
        met=directory("resources/tools/metilene"),
        met_out="results/dmr_calls/{group2}/metilene_output.bed",
    output:
        "results/dmr_calls/{group2}/plots/dmr_qval.0.05.bedgraph",
        report(
            "results/dmr_calls/{group2}/plots/dmr_qval.0.05.pdf",
            caption="../report/metilene_plots.rst",
            category="DMR plots",
            subcategory="Metilene plots",
            labels=lambda wildcards: {
                "experiment 1": config["ref_sample"],
                "experiment 2": wildcards.group2,
            },
        ),
    conda:
        "../envs/metilene.yaml"
    params:
        path_prefix=lambda wildcards: f"results/dmr_calls/{wildcards.group2}/plots/dmr",
        sample=lambda wildcards: {wildcards.group2.split("-")[0]},
        base_exp_number=config["ref_sample_number"],
    shell:
        """
        perl {input.met}/metilene_output.pl -q {input.met_out} -o {params.path_prefix} -a {params.base_exp_number} -b {wildcards.group2}
        """


# Still included for debugging reason but not used right now
rule _plot_pvals:
    input:
        met_out="results/dmr_calls/{group2}/metilene_output.bed",
    output:
        "results/dmr_calls/{group2}/plots/pvals.png",
    conda:
        "../envs/plot.yaml"
    params:
        path_prefix=lambda wildcards: f"results/dmr_calls/{wildcards.group2}",
        sample=lambda wildcards: {wildcards.group2.split("-")[0]},
    script:
        "../scripts/plot_pvals.py"
