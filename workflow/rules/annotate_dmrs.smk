rule download_regulatory_elements:
    output:
        "resources/ref/regulatory_elements.gff3",
    params:
        species=chromosome_conf["species"],
        species_cap=chromosome_conf["species"].capitalize(),
        build=chromosome_conf["build"],
        release=chromosome_conf["release"],
    shell:
        """
        wget -O {output}.gz https://ftp.ensembl.org/pub/release-{params.release}/regulation/{params.species}/{params.build}/annotation/{params.species_cap}.{params.build}.regulatory_features.v{params.release}.gff3.gz
        gzip -d {output}.gz
        """


rule annotate_regulatory_elements:
    input:
        metilene="results/{platform}/{caller}/dmr_calls/{group2}/metilene_output.bed",
        gene_annotation="resources/ref/regulatory_elements.gff3",
    output:
        "results/{platform}/{caller}/dmr_calls/{group2}/regulatory_elements/regulatory_elements.tsv",
    conda:
        "../envs/bedtools.yaml"
    shell:
        """
        bedtools intersect -a {input.metilene} -b {input.gene_annotation} -wa -wb > {output}
        """


rule add_regulatory_elements_header:
    input:
        "results/{platform}/{caller}/dmr_calls/{group2}/regulatory_elements/regulatory_elements.tsv",
    output:
        "results/{platform}/{caller}/dmr_calls/{group2}/regulatory_elements/regulatory_elements_complete.tsv",
    shell:
        """
        echo -e "chr\tstart_dmr\tend_dmr\tq-value\tmean_methylation_difference\tnumber_CpGs\tp(MWU)\tp(2DKS)\tmean_g1\tmean_g2\tseqif\tsource\ttype\tstart_feature\tend_feature\tscore\tstrand\tphase\tattributes" > {output}
        cat {input} >> {output}
        """


rule postprocess_regulatory_elements:
    input:
        "results/{platform}/{caller}/dmr_calls/{group2}/regulatory_elements/regulatory_elements_complete.tsv",
    output:
        "results/{platform}/{caller}/dmr_calls/{group2}/regulatory_elements/regulatory_elements_postprocessed.tsv",
    conda:
        "../envs/python_standard.yaml"
    script:
        "../scripts/postprocess_regulatory_elements.py"


rule get_gene_elements_annotation:
    output:
        "resources/ref/annotation.gtf.gz",
    params:
        species=chromosome_conf["species"],
        build=chromosome_conf["build"],
        release=chromosome_conf["release"],
    log:
        "logs/get_annotation.log",
    cache: "omit-software"  # save space and time with between workflow caching (see docs)
    wrapper:
        "v3.3.5/bio/reference/ensembl-annotation"


rule generate_txdb:
    input:
        "resources/ref/annotation.gtf.gz",
    output:
        txdb="resources/ref/txdb.db",
        txnames="resources/ref/txnames.tsv",
    log:
        "logs/generate_txdb.log",
    conda:
        "../envs/genomicfeatures.yaml"
    script:
        "../scripts/generate_txdb.R"


rule annotate_gene_elements:
    input:
        metilene="results/{platform}/{caller}/dmr_calls/{group2}/metilene_output.bed",
        txdb="resources/ref/txdb.db",
        txnames="resources/ref/txnames.tsv",
    output:
        chipseeker="results/{platform}/{caller}/dmr_calls/{group2}/genes_transcripts/chipseeker.tsv",
    log:
        "logs/annotate_gene_elements_{platform}_{caller}_{group2}.log",
    conda:
        "../envs/chipseeker.yaml"
    script:
        "../scripts/chipseeker.R"


rule rename_chipseeker_column:
    input:
        "data/input.tsv",
    output:
        "data/output.tsv",
    shell:
        """
        awk 'NR==1{{gsub("q_value", "qval")}}1' {input} > {output}
        """


rule datavzrd_annotations:
    input:
        config=workflow.source_path("../resources/dmrs_annotated.yaml"),
        genes_transcripts="results/{platform}/{caller}/dmr_calls/{group2}/genes_transcripts/chipseeker_postprocessed.tsv",
        regulatory_elements="results/{platform}/{caller}/dmr_calls/{group2}/regulatory_elements/regulatory_elements_postprocessed.tsv",
    output:
        report(
            directory("results/{platform}/{caller}/datavzrd-report/{group2}"),
            caption="../report/annotations.rst",
            htmlindex="index.html",
            category="Annotated DMRs",
            subcategory=lambda wildcards: f"{wildcards.platform} - {wildcards.caller}",
            labels=lambda wildcards: {
                "experiment 1": config["ref_sample"],
                "experiment 2": wildcards.group2,
            },
        ),
    params:
        base_experiment=lambda wildcards: config["ref_sample"],
    wrapper:
        "v3.13.0/utils/datavzrd"
