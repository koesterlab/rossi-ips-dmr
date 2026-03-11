rule get_old_rna_seq_fastqs:
    output:
        "rna_old/resources/fastqs/{accession}.fastq.gz",
    log:
        "logs/get_old_rna_seq_fastqs/{accession}.log",
    params:
        extra="--skip-technical",
    threads: 6
    wrapper:
        "v7.6.0/bio/sra-tools/fasterq-dump"


rule prepare_kallisto_sleuth_old_rna_seq:
    input:
        fastqs=expand(
            "rna_old/resources/fastqs/{accession}.fastq.gz",
            accession=config["rna_accessions_old"].keys(),
        ),
    output:
        samples="rna_old/config/samples_old.tsv",
        units="rna_old/config/units_old.tsv",
    conda:
        "../envs/python.yaml"
    log:
        "logs/prepare_kallisto_sleuth_old.log",
    params:
        types=config["rna_accessions_old"],
    script:
        "../scripts/prepare_kallisto_sleuth_old.py"


rule unzip_rna_new:
    input:
        "rna_new/KOLF_Trilineage_RNAseq.zip",
    output:
        expand(
            "rna_new/resources/fastqs/{barcode}",
            barcode=[
                f"SQK-NBD114-24_barcode{str(i).zfill(2)}.bam" for i in range(1, 13)
            ],
        ),
    log:
        "logs/unzip_rna_new.log",
    threads: 4
    shell:
        """
        UNZIP_DISABLE_ZIPBOMB_DETECTION=TRUE unzip -o {input} -d $(dirname {output[0]}) > {log} 2>&1
        """


rule rna_bam_to_fastq_rna_new:
    input:
        bam="rna_new/resources/fastqs/{sample}.bam",
    output:
        fq="rna_new/resources/fastqs/{sample}.fastq.gz",
    conda:
        "../envs/samtools.yaml"
    log:
        "logs/bam_to_fastq/{sample}.log",
    shell:
        """
        samtools fastq {input.bam} | gzip > {output.fq} 2> {log}
        """


rule prepare_kallisto_sleuth_new_rna_seq:
    input:
        fastqs=expand(
            "rna_new/resources/fastqs/{sample}.fastq.gz",
            sample=config["rna_accessions_new"].keys(),
        ),
    output:
        samples="rna_new/config/samples_new.tsv",
        units="rna_new/config/units_new.tsv",
    conda:
        "../envs/python.yaml"
    log:
        "logs/prepare_kallisto_sleuth_new.log",
    params:
        types=config["rna_accessions_new"],
    script:
        "../scripts/prepare_kallisto_sleuth_new.py"


# ---------------------------------------------------------------------------
# Compare differential methylation (DMR) to differential expression
# ---------------------------------------------------------------------------
# The {rna_data} wildcard is always "rna_old" or "rna_new".
# The kallisto-sleuth modules are prefixed with rna_{rna_data}/base_{base},
# so the diffexp tables live at:
#   rna_old/base_{base}/results/tables/diffexp/condition.genes-representative.diffexp_postprocessed.tsv
#   rna_new/base_{base}/results/tables/diffexp/condition.genes-representative.diffexp_postprocessed.tsv


rule compare_dmr_to_diffexp_no_tfs:
    input:
        diffexp=lambda wildcards: (
            f"{wildcards.rna_data}/base_{wildcards.base}"
            "/results/tables/diffexp/"
            "condition.genes-representative.diffexp_postprocessed.tsv"
        ),
        # For a given base, ChIPseeker annotations exist for the 3 non-base layers.
        layer1=lambda wildcards: (
            "results/{platform}/{caller}/base_{base}/dmr_calls/"
            + get_non_base_layers(wildcards.base)[0]
            + "/genes_transcripts/chipseeker_postprocessed.tsv"
        ).format(**wildcards),
        layer2=lambda wildcards: (
            "results/{platform}/{caller}/base_{base}/dmr_calls/"
            + get_non_base_layers(wildcards.base)[1]
            + "/genes_transcripts/chipseeker_postprocessed.tsv"
        ).format(**wildcards),
        layer3=lambda wildcards: (
            "results/{platform}/{caller}/base_{base}/dmr_calls/"
            + get_non_base_layers(wildcards.base)[2]
            + "/genes_transcripts/chipseeker_postprocessed.tsv"
        ).format(**wildcards),
    output:
        tsv="results/{platform}/{caller}/base_{base}/{rna_data}/diffexp_vs_dmrs_{annotation_type}.tsv",
        html="results/{platform}/{caller}/base_{base}/{rna_data}/diffexp_vs_dmrs_{annotation_type}.html",
    wildcard_constraints:
        # rna_data is strictly rna_old or rna_new – no slashes
        rna_data="rna_old|rna_new",
    conda:
        "../envs/python.yaml"
    params:
        annotation_type=lambda wildcards: wildcards.annotation_type,
        base=lambda wildcards: wildcards.base,
        # Pass the three non-base layer names so the script can label them
        # and look up the correct diffexp columns dynamically.
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
    wildcard_constraints:
        rna_data="rna_old|rna_new",
    conda:
        "../envs/python.yaml"
    params:
        annotation_type=lambda wildcards: wildcards.annotation_type,
    log:
        "logs/compare_dmr_to_diffexp_with_tfs/{platform}_{caller}_{base}_{rna_data}_{annotation_type}_tfs.log",
    script:
        "../scripts/compare_dmr_to_diffexp_with_tfs.py"


rule datavzrd_dmr_vs_diffexp_no_tfs:
    input:
        config=workflow.source_path("../resources/dmr_vs_diffexp_no_tfs.yaml"),
        table="results/{platform}/{caller}/base_{base}/{rna_data}/diffexp_vs_dmrs_{annotation_type}.tsv",
    output:
        report(
            directory(
                "results/{platform}/{caller}/base_{base}/{rna_data}/diffexp_vs_dmrs_no_tfs_{annotation_type}"
            ),
            caption="../report/diffexp_vs_dmrs.rst",
            htmlindex="index.html",
            category="DiffExp-DMRs Comparison",
            subcategory=lambda wildcards: "Comparisons",
            labels=lambda wildcards: {
                "base": wildcards.base,
                "rna_data": wildcards.rna_data,
                "type": "no transcription factors",
            },
        ),
    wildcard_constraints:
        rna_data="rna_old|rna_new",
    log:
        "logs/diffexp_dmvzrd/diffexp_dmr_datavzrd/{platform}_{caller}_{base}_{rna_data}_{annotation_type}.log",
    conda:
        "../envs/datavzrd.yaml"
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
            category="DiffExp-DMRs Comparison",
            subcategory=lambda wildcards: "Comparisons",
            labels=lambda wildcards: {
                "base": wildcards.base,
                "rna_data": wildcards.rna_data,
                "type": "with transcription factors",
            },
        ),
    wildcard_constraints:
        rna_data="rna_old|rna_new",
    log:
        "logs/diffexp_dmvzrd/diffexp_dmr_datavzrd_with_tfs/{platform}_{caller}_{base}_{rna_data}_{annotation_type}.log",
    wrapper:
        "v9.2.0/utils/datavzrd"
