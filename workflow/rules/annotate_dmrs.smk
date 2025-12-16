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
        metilene="results/{platform}/{caller}/dmr_calls/{group2}/metilene_output_focused.bed",
        gene_annotation="resources/ref/regulatory_elements.gff3",
    output:
        "results/{platform}/{caller}/dmr_calls/{group2}/regulatory_elements/regulatory_elements.tsv",
    conda:
        "../envs/bedtools.yaml"
    log:
        "logs/annotate_regulatory_elements/{platform}_{caller}_{group2}.log",
    shell:
        """
        bedtools intersect -a {input.metilene} -b {input.gene_annotation} -wa -wb > {output} 2> {log}
        """


rule add_regulatory_elements_header:
    input:
        "results/{platform}/{caller}/dmr_calls/{group2}/regulatory_elements/regulatory_elements.tsv",
    output:
        "results/{platform}/{caller}/dmr_calls/{group2}/regulatory_elements/regulatory_elements_complete.tsv",
    log:
        "logs/add_regulatory_elements_header/{platform}_{caller}_{group2}.log",
    shell:
        """
        echo -e "chr\tstart_dmr\tend_dmr\tq-value\tmean_methylation_difference\tnumber_CpGs\tp(MWU)\tp(2DKS)\tmean_g1\tmean_g2\tseqif\tsource\ttype\tstart_feature\tend_feature\tscore\tstrand\tphase\tattributes" > {output} 2> {log}
        cat {input} >> {output} 2> {log}
        """

# Prepare table for better readability with datavzrd
rule postprocess_regulatory_elements:
    input:
        "results/{platform}/{caller}/dmr_calls/{group2}/regulatory_elements/regulatory_elements_complete.tsv",
    output:
        "results/{platform}/{caller}/dmr_calls/{group2}/regulatory_elements/regulatory_elements_postprocessed.tsv",
    conda:
        "../envs/python_standard.yaml"
    log:
        "logs/postprocess_regulatory_elements/{platform}_{caller}_{group2}.log",
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
        "logs/get_gene_elements_annotation.log",
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


rule annotate_gene_elements_filtered:
    input:
        metilene="results/{platform}/{caller}/dmr_calls/{group2}/metilene_output_focused.bed",
        txdb="resources/ref/txdb.db",
        txnames="resources/ref/txnames.tsv",
    output:
        chipseeker="results/{platform}/{caller}/dmr_calls/{group2}/genes_transcripts/chipseeker_filtered.tsv",
    log:
        "logs/annotate_gene_elements_filtered/{platform}_{caller}_{group2}.log",
    conda:
        "../envs/chipseeker.yaml"
    script:
        "../scripts/chipseeker.R"


# We need the complete unfiltered list for rna seq comparison
rule annotate_gene_elements_complete:
    input:
        metilene="results/{platform}/{caller}/dmr_calls/{group2}/metilene_output.bed",
        txdb="resources/ref/txdb.db",
        txnames="resources/ref/txnames.tsv",
    output:
        chipseeker="results/{platform}/{caller}/dmr_calls/{group2}/genes_transcripts/chipseeker_complete.tsv",
    log:
        "logs/annotate_gene_elements_complete/{platform}_{caller}_{group2}.log",
    conda:
        "../envs/chipseeker.yaml"
    script:
        "../scripts/chipseeker.R"



rule get_ensembl_genes:
    input:
        "results/{platform}/{caller}/dmr_calls/{group2}/genes_transcripts/chipseeker_complete.tsv",
    output:
        "results/{platform}/{caller}/dmr_calls/{group2}/genes_transcripts/ensembl_genes.tsv",
    conda:
        "../envs/biomart.yaml"
    params:
        species=get_bioc_species_name(),
        version=config["resources"]["ref"]["release"],
    log:
        "logs/get_ensembl_genes/{platform}_{caller}_{group2}.log",
    script:
        "../scripts/get_ensembl_genes.R"



# We only merge with polars because the join in get_ensembl_genes.R does not work
rule annotate_chipseeker:
    input:
        chipseeker="results/{platform}/{caller}/dmr_calls/{group2}/genes_transcripts/chipseeker_{type}.tsv",
        genes="results/{platform}/{caller}/dmr_calls/{group2}/genes_transcripts/ensembl_genes.tsv",
    output:
        "results/{platform}/{caller}/dmr_calls/{group2}/genes_transcripts/chipseeker_postprocessed_{type}.tsv",
    conda:
        "../envs/python_standard.yaml"
    log:
        "logs/annotate_chipseeker/{platform}_{caller}_{group2}_{type}.log",
    script:
        "../scripts/annotate_chipseeker.py"


rule datavzrd_annotations:
    input:
        config=workflow.source_path("../resources/dmrs_annotated.yaml"),
        genes_transcripts="results/{platform}/{caller}/dmr_calls/{group2}/genes_transcripts/chipseeker_postprocessed_filtered.tsv",
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
    log:
        "logs/datavzrd_annotations/{platform}_{caller}_{group2}.log",
    wrapper:
        "v3.13.0/utils/datavzrd"
