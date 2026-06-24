############################################ Annotate regulatory elements ############################################


rule download_regulatory_elements:
    output:
        "resources/ref/regulatory_elements.gff3",
    params:
        species=chromosome_conf["species"],
        species_cap=chromosome_conf["species"].capitalize(),
        build=chromosome_conf["build"],
        release=chromosome_conf["release"],
    log:
        "logs/download_regulatory_elements.log",
    shell:
        """
        wget -O {output}.gz https://ftp.ensembl.org/pub/release-{params.release}/regulation/{params.species}/{params.build}/annotation/{params.species_cap}.{params.build}.regulatory_features.v{params.release}.gff3.gz 2> {log}
        gzip -d {output}.gz
        """


rule annotate_regulatory_elements:
    input:
        metilene="results/{platform}/{caller}/base_{base}/dmr_calls/{group2}/metilene_output_focused_0.05.bed",
        gene_annotation="resources/ref/regulatory_elements.gff3",
    output:
        "results/{platform}/{caller}/base_{base}/dmr_calls/{group2}/regulatory_elements/regulatory_elements.tsv",
    conda:
        "../envs/bedtools.yaml"
    log:
        "logs/annotate_regulatory_elements/{platform}_{caller}_{base}_{group2}.log",
    shell:
        """
        bedtools intersect -a {input.metilene} -b {input.gene_annotation} -wa -wb > {output} 2> {log}
        """


rule add_regulatory_elements_header:
    input:
        "results/{platform}/{caller}/base_{base}/dmr_calls/{group2}/regulatory_elements/regulatory_elements.tsv",
    output:
        "results/{platform}/{caller}/base_{base}/dmr_calls/{group2}/regulatory_elements/regulatory_elements_complete.tsv",
    log:
        "logs/add_regulatory_elements_header/{platform}_{caller}_{base}_{group2}.log",
    shell:
        """
        echo -e "chr\tstart_dmr\tend_dmr\tq-value\tmean_methylation_difference\tnumber_CpGs\tp(MWU)\tp(2DKS)\tmean_g1\tmean_g2\tseqif\tsource\ttype\tstart_feature\tend_feature\tscore\tstrand\tphase\tattributes" > {output} 2> {log}
        cat {input} >> {output} 2> {log}
        """


# Prepare table for better readability with datavzrd
rule postprocess_regulatory_elements:
    input:
        "results/{platform}/{caller}/base_{base}/dmr_calls/{group2}/regulatory_elements/regulatory_elements_complete.tsv",
    output:
        "results/{platform}/{caller}/base_{base}/dmr_calls/{group2}/regulatory_elements/regulatory_elements_postprocessed.tsv",
    conda:
        "../envs/python.yaml"
    log:
        "logs/postprocess_regulatory_elements/{platform}_{caller}_{base}_{group2}.log",
    script:
        "../scripts/postprocess_regulatory_elements.py"


############################################### Annotate gene elements ############################################


# Gene elements are e.g. promoters, exons, introns, etc.
rule get_gene_elements_annotation:
    output:
        "resources/ref/annotation.gtf.gz",
    params:
        species=chromosome_conf["species"],
        build=chromosome_conf["build"],
        release=chromosome_conf["release"],
    log:
        "logs/get_gene_elements_annotation.log",
    cache: "omit-software"  # save space and time with between workflow caching (see docs)
    wrapper:
        "v3.3.5/bio/reference/ensembl-annotation"


rule generate_txdb_from_gene_elements:
    input:
        "resources/ref/annotation.gtf.gz",
    output:
        txdb="resources/ref/txdb.db",
        txnames="resources/ref/txnames.rds",
    log:
        "logs/generate_txdb_from_gene_elements.log",
    conda:
        "../envs/genomicfeatures.yaml"
    script:
        "../scripts/generate_txdb.R"


rule annotate_dmrs_with_gene_elements:
    input:
        metilene="results/{platform}/{caller}/base_{base}/dmr_calls/{group2}/metilene_output_focused_{fdr}.bed",
        txdb="resources/ref/txdb.db",
        txnames="resources/ref/txnames.rds",
    output:
        chipseeker="results/{platform}/{caller}/base_{base}/dmr_calls/{group2}/genes_transcripts/chipseeker_{fdr}.tsv",
    log:
        "logs/annotate_dmrs_with_gene_elements/{platform}_{caller}_{base}_{group2}.log",
    conda:
        "../envs/chipseeker.yaml"
    script:
        "../scripts/chipseeker.R"


# We want real gene names like SOX2 instead of Ensembl transcript IDs.
rule get_ensembl_gene_names_from_dmrs:
    input:
        "results/{platform}/{caller}/base_{base}/dmr_calls/{group2}/genes_transcripts/chipseeker_{fdr}.tsv",
    output:
        "results/{platform}/{caller}/base_{base}/dmr_calls/{group2}/genes_transcripts/ensembl_genes_{fdr}.tsv",
    conda:
        "../envs/biomart.yaml"
    params:
        species=get_bioc_species_name(),
        version=config["resources"]["ref"]["release"],
    log:
        "logs/get_ensembl_gene_names_from_dmrs/{platform}_{caller}_{base}_{group2}.log",
    # Use unrealistc high memory to avoid parallel computation since ensembl then detects DOS attacks
    resources:
        mem_mb=16000,
    script:
        "../scripts/get_ensembl_genes.R"


rule annotate_dmrs_with_ensembl_gene_names:
    input:
        chipseeker="results/{platform}/{caller}/base_{base}/dmr_calls/{group2}/genes_transcripts/chipseeker_{fdr}.tsv",
        genes="results/{platform}/{caller}/base_{base}/dmr_calls/{group2}/genes_transcripts/ensembl_genes_{fdr}.tsv",
    output:
        "results/{platform}/{caller}/base_{base}/dmr_calls/{group2}/genes_transcripts/chipseeker_postprocessed_{fdr}.tsv",
    conda:
        "../envs/python.yaml"
    log:
        "logs/annotate_dmrs_with_ensembl_gene_names/{platform}_{caller}_{base}_{group2}.log",
    script:
        "../scripts/annotate_chipseeker.py"


rule dmr_heatmap:
    input:
        lambda wildcards: expand(
            "results/{{platform}}/{{caller}}/base_{{base}}/dmr_calls/{group2}/genes_transcripts/chipseeker_postprocessed_0.05.tsv",
            group2=get_non_base_layers(wildcards.base),
        ),
    output:
        report(
            "results/{platform}/{caller}/base_{base}/dmr_calls/heatmaps/{type}.png",
            caption="../report/heatmap.rst",
            category="DMR plots",
            subcategory=lambda wildcards: f"Heatmaps: {wildcards.platform}",
            labels=lambda wildcards: {
                "base": wildcards.base,
                "genetic element": wildcards.type,
            },
        ),
    conda:
        "../envs/plot.yaml"
    log:
        "logs/dmr_heatmap/{platform}_{caller}_{base}_{type}.log",
    resources:
        mem_mb=16000,
    params:
        base = lambda wildcards: wildcards.base
    script:
        "../scripts/dmr-heatmap.py"


rule datavzrd_annotations:
    input:
        config=workflow.source_path("../resources/dmrs_annotated.yaml"),
        genes_transcripts="results/{platform}/{caller}/base_{base}/dmr_calls/{group2}/genes_transcripts/chipseeker_postprocessed_0.05.tsv",
        regulatory_elements="results/{platform}/{caller}/base_{base}/dmr_calls/{group2}/regulatory_elements/regulatory_elements_postprocessed.tsv",
    output:
        report(
            directory("results/{platform}/{caller}/base_{base}/dmr_calls/datavzrd-report/{group2}"),
            caption="../report/annotations.rst",
            htmlindex="index.html",
            category="Annotated DMRs",
            subcategory=lambda wildcards: f"{wildcards.platform} - Base: {wildcards.base}",
            labels=lambda wildcards: {
                "comparison": wildcards.group2,
            },
        ),
    params:
        base_experiment=lambda wildcards: wildcards.base,
    log:
        "logs/datavzrd_annotations/{platform}_{caller}_{base}_{group2}.log",
    wrapper:
        "641c90c4da86d4acf2022f347f3c8017334c0f44/utils/datavzrd"
