
rule download_regulations:
    output:
        "resources/ref/regulations.gff3",
    shell:
        """
        wget -O {output}.gz https://ftp.ensembl.org/pub/release-112/regulation/homo_sapiens/GRCh38/annotation/Homo_sapiens.GRCh38.regulatory_features.v112.gff3.gz
        gzip -d {output}.gz
        """


rule annotate_regulations_dmrs:
    input:
        metilene="results/dmr_calls/{group2}/metilene_output.bed",
        gene_annotation="resources/ref/regulations.gff3",
    output:
        "results/dmr_calls/{group2}/reculations/annotated.tsv",
    conda:
        "../envs/bedtools.yaml"
    shell:
        """
        bedtools intersect -a {input.metilene} -b {input.gene_annotation} -wa -wb > {output}
        """


rule add_regulations_header:
    input:
        "results/dmr_calls/{group2}/reculations/annotated.tsv",
    output:
        "results/dmr_calls/{group2}/reculations/annotated_header.tsv",
    shell:
        """
        echo -e "chr\tstart_dmr\tend_dmr\tq-value\tmean_methylation_difference\tnumber_CpGs\tp(MWU)\tp(2DKS)\tmean_g1\tmean_g2\tseqif\tsource\ttype\tstart_feature\tend_feature\tscore\tstrand\tphase\tattributes" > {output}
        cat {input} >> {output}
        """


rule get_annotation:
    output:
        "resources/ref/annotation.gtf.gz",
    params:
        species=config_ref["species"],
        build=config_ref["build"],
        release=config_ref["release"],
    log:
        "logs/get_annotation.log",
    cache: "omit-software"  # save space and time with between workflow caching (see docs)
    wrapper:
        "v3.3.5/bio/reference/ensembl-annotation"


rule generate_txdb:
    input:
        "resources/ref/annotation.gtf.gz",
    output:
        txdb="resources/txdb.db",
        txnames="resources/txnames.tsv",
    log:
        "logs/generate_txdb.log",
    conda:
        "../envs/genomicfeatures.yaml"
    script:
        "../scripts/generate_txdb.R"


rule chipseeker_annotate:
    input:
        metilene="results/dmr_calls/{group2}/metilene_output.bed",
        txdb="resources/txdb.db",
        txnames="resources/txnames.tsv",
    output:
        chipseeker="results/dmr_calls/{group2}/genes_transcripts/chipseeker_regions.tsv",
    log:
        "logs/chipseeker_annotate_{group2}.log",
    conda:
        "../envs/chipseeker.yaml"
    script:
        "../scripts/chipseeker.R"


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
        config=workflow.source_path("../resources/datavzrd.yaml"),
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
