config_ref = config["resources"]["ref"]


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
        chipseeker="results/dmr_calls/{group2}/chipseeker_regions.tsv",
    log:
        "logs/chipseeker_annotate_{group2}.log",
    conda:
        "../envs/chipseeker.yaml"
    script:
        "../scripts/chipseeker.R"


rule create_heatmap:
    input:
        expand(
            "results/dmr_calls/{sample}-ref-sorted/chipseeker_regions.tsv",
            sample=["BC02", "BC03", "BC04"],
        ),
    output:
        "results/dmr_calls/heatmaps/{type}.png",
    conda:
        "../envs/plot.yaml"
    script:
        "../scripts/heatmap.py"
