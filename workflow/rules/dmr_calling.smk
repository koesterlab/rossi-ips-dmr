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
        varlo_path=config["varlo_path"],
        pipeline_path=config["pipeline_path"],
        event="PRESENT",
    shell:
        """
        cd {params.varlo_path}
        cargo run --release -- filter-calls control-fdr --mode local-smart {params.pipeline_path}{input} --events {params.event} --fdr 0.01 > {params.pipeline_path}{output}
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
            category="DMR plots",
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


rule plot_pvals:
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


rule download_genes_and_transcripts:
    output:
        "resources/genes_transcripts.gff3",
    shell:
        """
        wget -O {output}.gz https://ftp.ensembl.org/pub/release-112/gff3/homo_sapiens/Homo_sapiens.GRCh38.112.gff3.gz &&
        gzip -d {output}.gz
        """


rule download_regulations:
    output:
        "resources/regulations.gff3",
    shell:
        """
        wget -O {output}.gz https://ftp.ensembl.org/pub/release-112/regulation/homo_sapiens/GRCh38/annotation/Homo_sapiens.GRCh38.regulatory_features.v112.gff3.gz
        gzip -d {output}.gz
        """


rule annotate_dmrs:
    input:
        metilene="results/dmr_calls/{group2}/metilene_output.bed",
        gene_annotation="resources/{comparison}.gff3",
    output:
        "results/dmr_calls/{group2}/{comparison}/annotated.tsv",
    conda:
        "../envs/bedtools.yaml"
    shell:
        """
        bedtools intersect -a {input.metilene} -b {input.gene_annotation} -wa -wb > {output}
        """


rule add_header:
    input:
        "results/dmr_calls/{group2}/{comparison}/annotated.tsv",
    output:
        "results/dmr_calls/{group2}/{comparison}/annotated_header.tsv",
    shell:
        """
        echo -e "chr\tstart_dmr\tend_dmr\tq-value\tmean_methylation_difference\tnumber_CpGs\tp(MWU)\tp(2DKS)\tmean_g1\tmean_g2\tseqif\tsource\ttype\tstart_feature\tend_feature\tscore\tstrand\tphase\tattributes" > {output}
        cat {input} >> {output}
        """


rule postprocess_annotation:
    input:
        "results/dmr_calls/{group2}/{comparison}/annotated_header.tsv",
    output:
        "results/dmr_calls/{group2}/{comparison}/annotated_complete.tsv",
    conda:
        "../envs/pandas.yaml"
    script:
        "../scripts/postprocess_annotation.py"


rule datavzrd:
    input:
        config="resources/datavzrd.yaml",
        genes_transcripts="results/dmr_calls/{group2}/genes_transcripts/annotated_complete.tsv",
        regulations="results/dmr_calls/{group2}/regulations/annotated_complete.tsv",
    output:
        report(
            directory("results/datavzrd-report/{group2}"),
            htmlindex="index.html",
            category="Annotated DMRs",
            labels=lambda wildcards: {
                "experiment 1": config["ref_sample"],
                "experiment 2": wildcards.group2,
            },
        ),
    params:
        base_experiment=lambda wildcards: config["ref_sample"],
    wrapper:
        "v3.13.0/utils/datavzrd"
