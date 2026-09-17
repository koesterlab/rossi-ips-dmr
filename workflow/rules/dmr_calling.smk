rule download_metilene:
    output:
        temp("resources/tools/metilene_v02-9.tar.gz"),
    log:
        "logs/download_metilene.log",
    shell:
        "wget -O {output} http://www.bioinf.uni-leipzig.de/Software/metilene/metilene_v02-9.tar.gz 2> {log}"


rule unpack_metilene:
    input:
        "resources/tools/metilene_v02-9.tar.gz",
    output:
        directory("resources/tools/metilene"),
    log:
        "logs/unpack_metilene.log",
    shell:
        """
        mkdir -p {output}
        tar -xzf {input} -C {output} --strip-components=1 2> {log}
        """


rule metilene_input:
    input:
        "results/{platform}/{caller}/meth_calling/calls_{fdr}.parquet",
    output:
        "results/{platform}/{caller}/base_{base}/dmr_calls/{group2}/{fdr}/metilene_input.txt",
    conda:
        "../envs/python.yaml"
    resources:
        mem_mb=32000,
    log:
        "logs/metilene_input/{platform}_{caller}_{base}_{group2}_{fdr}.log",
    script:
        "../scripts/metilene_input.py"

# Output cols:
# | chr | start | stop | q-value | mean methylation difference | #CpGs | p (MWU) | p (2D KS) | mean g1 | mean g2 |
rule call_metilene:
    input:
        "results/{platform}/{caller}/base_{base}/dmr_calls/{group2}/{fdr}/metilene_input.txt",
    output:
        "results/{platform}/{caller}/base_{base}/dmr_calls/{group2}/{fdr}/metilene_output.bed",
    conda:
        "../envs/metilene.yaml"
    threads: 4
    log:
        "logs/call_metilene/{platform}_{caller}_{base}_{group2}_{fdr}.log",
    shell:
        "metilene -d 0.01 -t {threads} -m 10 -a {wildcards.group2} -b {wildcards.base} {input} > {output} 2> {log}"


rule focus_dmrs:
    """Keep only DMRs of {germ_layer} that do not overlap DMRs of the other three non-base layers."""
    input:
        this="results/{platform}/{caller}/base_{base}/dmr_calls/{germ_layer}/{fdr}/metilene_output.bed",
        others=lambda wildcards: expand(
            "results/{{platform}}/{{caller}}/base_{{base}}/dmr_calls/{layer}/{{fdr}}/metilene_output.bed",
            layer=[
                layer
                for layer in get_non_base_layers(wildcards.base)
                if layer != wildcards.germ_layer
            ],
        ),
    output:
        "results/{platform}/{caller}/base_{base}/dmr_calls/{germ_layer}/{fdr}/metilene_output_focused.bed",
    conda:
        "../envs/bedtools.yaml"
    log:
        "logs/focus_dmrs/{platform}_{caller}_base_{base}_{germ_layer}_{fdr}.log",
    script:
        "../scripts/focus_dmrs.py"


rule metilene_plots:
    input:
        met="resources/tools/metilene",
        met_out="results/{platform}/{caller}/base_{base}/dmr_calls/{group2}/0.05/metilene_output_focused.bed",
    output:
        bed="results/{platform}/{caller}/base_{base}/dmr_calls/{group2}/plots/dmr_qval.0.05.bedgraph",
        pdf=report(
            "results/{platform}/{caller}/base_{base}/dmr_calls/{group2}/plots/dmr_qval.0.05.pdf",
            caption="../report/metilene_plots.rst",
            category="DMR plots",
            subcategory=lambda wildcards: f"Metilene: {wildcards.platform} - {wildcards.caller}",
            labels=lambda wildcards: {
                "base": wildcards.base,
                "layer": wildcards.group2,
            },
        ),
    conda:
        "../envs/metilene.yaml"
    log:
        "logs/metilene_plots/{platform}_{caller}_{base}_{group2}.log",
    shell:
        # Group order must match the metilene call above (-a group2 -b base).
        "perl {input.met}/metilene_output.pl -q {input.met_out} -o $(dirname {output.bed})/dmr -a {wildcards.group2} -b {wildcards.base} 2> {log}"
