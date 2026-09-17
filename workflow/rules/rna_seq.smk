rule get_old_rna_seq_fastqs:
    output:
        "resources/rna_seq_old/fastqs/{accession}.fastq.gz",
    log:
        "logs/get_old_rna_seq_fastqs/{accession}.log",
    params:
        extra="--skip-technical",
    threads: 6
    wrapper:
        "v9.4.0/bio/sra-tools/fasterq-dump"

rule unzip_rna_new:
    input:
        "resources/rna_seq_new/KOLF_Trilineage_RNAseq_new.zip",
    output:
        expand(
            "resources/rna_seq_new/bams/{barcode}.bam",
            barcode=config["rna_accessions_new"],
        ),
    log:
        "logs/unzip_rna_new.log",
    threads: 4
    shell:
        """
        mkdir -p $(dirname {output[0]})
        UNZIP_DISABLE_ZIPBOMB_DETECTION=TRUE unzip -o {input} -d $(dirname {output[0]}) > {log} 2>&1
        """


rule rna_bam_to_fastq_rna_new:
    input:
        bam="resources/rna_seq_new/bams/{barcode}.bam",
    output:
        fq="resources/rna_seq_new/fastqs/{barcode}.fastq.gz",
    conda:
        "../envs/samtools.yaml"
    log:
        "logs/bam_to_fastq/{barcode}.log",
    shell:
        "samtools fastq {input.bam} 2> {log} | gzip > {output.fq}"


rule prepare_kallisto_sleuth:
    input:
        old_fastqs=expand(
            "resources/rna_seq_old/fastqs/{accession}.fastq.gz",
            accession=config["rna_accessions_old"],
        ),
        new_fastqs=expand(
            "resources/rna_seq_new/fastqs/{accession}.fastq.gz",
            accession=config["rna_accessions_new"],
        ),
    output:
        samples="config/samples.tsv",
        units="config/units.tsv",
    conda:
        "../envs/python.yaml"
    params:
        labels_old=config["rna_accessions_old"],
        labels_new=config["rna_accessions_new"],
    log:
        "logs/prepare_kallisto_sleuth.log",
    script:
        "../scripts/prepare_kallisto_sleuth.py"


# ---------------------------------------------------------------------------
# Compare differential methylation (DMRs) to differential expression
# ---------------------------------------------------------------------------


def short_layer(layer):
    """Germ layer name as used in the diffexp model names, e.g. "ectoderm" -> "ecto"."""
    return layer.removesuffix("derm")


def diffexp_models(wildcards):
    """
    For each non-base layer, return the configured sleuth model comparing it to
    the base, together with the sign needed to orient its effect sizes as
    "layer vs. base". Only one direction of each contrast is defined in the
    config; if it is "base vs. layer", the effect sizes have to be flipped.
    """
    models = config["kallisto_sleuth"]["diffexp"]["models"]
    dataset = wildcards.rna_data.removeprefix("rna_")
    base = short_layer(wildcards.base)
    result = []
    for layer in map(short_layer, get_non_base_layers(wildcards.base)):
        model = f"{layer}_vs_{base}_{dataset}"
        if model in models:
            result.append((model, 1))
        else:
            result.append((f"{base}_vs_{layer}_{dataset}", -1))
    return result


def diffexp_tables(wildcards):
    return [
        f"results/tables/diffexp/{model}.genes-representative.diffexp_postprocessed.tsv"
        for model, _ in diffexp_models(wildcards)
    ]


def diffexp_signs(wildcards):
    return [sign for _, sign in diffexp_models(wildcards)]


rule compare_dmr_to_diffexp_no_tfs:
    input:
        # Both lists are ordered like get_non_base_layers(base)
        diffexp=diffexp_tables,
        dmrs=chipseeker_tables,
        val_genes="resources/rna_seq/val_genes.tsv",
    output:
        tsv="results/{platform}/{caller}/base_{base}/{rna_data}/diffexp_vs_dmrs_{annotation_type}.tsv",
        dmr_diffexp="results/{platform}/{caller}/base_{base}/{rna_data}/diffexp_vs_dmrs_{annotation_type}.pdf",
        val_genes="results/{platform}/{caller}/base_{base}/{rna_data}/val_genes_{annotation_type}.tsv",
    conda:
        "../envs/python.yaml"
    params:
        diffexp_signs=diffexp_signs,
        non_base_layers=lambda wildcards: get_non_base_layers(wildcards.base),
    log:
        "logs/compare_dmr_to_diffexp_no_tfs/{platform}_{caller}_{base}_{rna_data}_{annotation_type}.log",
    script:
        "../scripts/compare_dmr_to_diffexp_no_tfs.py"


rule get_tfs_from_collectri:
    output:
        "resources/rna_seq/tf_list/collectri_tf_list.tsv",
    conda:
        "../envs/decoupler.yaml"
    log:
        "logs/get_tfs_from_collectri.log",
    script:
        "../scripts/get_tfs_from_collectri.py"


rule compare_dmr_to_diffexp_with_tfs:
    input:
        tf_list="resources/rna_seq/tf_list/collectri_tf_list.tsv",
        comp="results/{platform}/{caller}/base_{base}/{rna_data}/diffexp_vs_dmrs_{annotation_type}.tsv",
    output:
        focus_tfs="results/{platform}/{caller}/base_{base}/{rna_data}/tfs/diffexp_vs_dmrs_tfs_only_{annotation_type}.tsv",
        comp_tf_adj="results/{platform}/{caller}/base_{base}/{rna_data}/tfs/diffexp_vs_dmrs_{annotation_type}.tsv",
        plot="results/{platform}/{caller}/base_{base}/{rna_data}/tfs/diffexp_vs_dmrs_{annotation_type}.html",
    conda:
        "../envs/python.yaml"
    log:
        "logs/compare_dmr_to_diffexp_with_tfs/{platform}_{caller}_{base}_{rna_data}_{annotation_type}_tfs.log",
    script:
        "../scripts/compare_dmr_to_diffexp_with_tfs.py"


rule datavzrd_dmr_vs_diffexp_no_tfs:
    input:
        config=workflow.source_path("../resources/dmr_vs_diffexp_no_tfs.yaml"),
        table="results/{platform}/{caller}/base_{base}/{rna_data}/diffexp_vs_dmrs_{annotation_type}.tsv",
        val_genes="results/{platform}/{caller}/base_{base}/{rna_data}/val_genes_{annotation_type}.tsv",
    output:
        report(
            directory(
                "results/{platform}/{caller}/base_{base}/{rna_data}/diffexp_vs_dmrs_no_tfs_{annotation_type}"
            ),
            caption="../report/diffexp_vs_dmrs.rst",
            htmlindex="index.html",
            category=lambda wildcards: f"DiffExp-DMRs Comparison - {wildcards.rna_data}",
            subcategory=lambda wildcards: f"{wildcards.platform} - No tfs",
            labels=lambda wildcards: {
                "base": wildcards.base,
            },
        ),
    log:
        "logs/diffexp_dmvzrd/diffexp_dmr_datavzrd/{platform}_{caller}_{base}_{rna_data}_{annotation_type}.log",
    wrapper:
        "v9.2.0/utils/datavzrd"


rule datavzrd_dmr_vs_diffexp_with_tfs:
    input:
        config=workflow.source_path("../resources/dmr_vs_diffexp_with_tfs.yaml"),
        complete="results/{platform}/{caller}/base_{base}/{rna_data}/tfs/diffexp_vs_dmrs_{annotation_type}.tsv",
        focus_tfs="results/{platform}/{caller}/base_{base}/{rna_data}/tfs/diffexp_vs_dmrs_tfs_only_{annotation_type}.tsv",
    output:
        report(
            directory(
                "results/{platform}/{caller}/base_{base}/{rna_data}/diffexp_vs_dmrs_with_tfs_{annotation_type}"
            ),
            caption="../report/diffexp_vs_dmrs.rst",
            htmlindex="index.html",
            category=lambda wildcards: f"DiffExp-DMRs Comparison - {wildcards.rna_data}",
            subcategory=lambda wildcards: f"{wildcards.platform} - With tfs",
            labels=lambda wildcards: {
                "base": wildcards.base,
            },
        ),
    log:
        "logs/diffexp_dmvzrd/diffexp_dmr_datavzrd_with_tfs/{platform}_{caller}_{base}_{rna_data}_{annotation_type}.log",
    wrapper:
        "v9.2.0/utils/datavzrd"


rule view_val_genes:
    input:
        config=workflow.source_path("../resources/val_genes.rds"),
    output:
        tsv="resources/rna_seq/val_genes.tsv",
    conda:
        "../envs/human_annotation.yaml"
    log:
        "logs/val_genes/val_genes.log",
    script:
        "../scripts/view_val_genes.R"
