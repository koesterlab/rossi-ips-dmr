rule get_old_rna_seq_fastqs:
    output:
        "resources/rna_seq_old/fastqs/{accession}.fastq.gz",
    log:
        "logs/get_old_rna_seq_fastqs/{accession}.log",
    params:
        extra="--skip-technical",
    threads: 6
    wrapper:
        "v7.6.0/bio/sra-tools/fasterq-dump"

rule unzip_rna_new:
    input:
        "resources/rna_seq_new/KOLF_Trilineage_RNAseq.zip",
    output:
        expand(
            "resources/rna_seq_new/bams/{barcode}",
            barcode=[
                f"SQK-NBD114-24_barcode{str(i).zfill(2)}.bam" for i in range(1, 13)
            ],
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
        bam="resources/rna_seq_new/bams/{sample}.bam",
    output:
        fq="resources/rna_seq_new/fastqs/{sample}.fastq.gz",
    conda:
        "../envs/samtools.yaml"
    log:
        "logs/bam_to_fastq/{sample}.log",
    shell:
        """
        mkdir -p $(dirname {output.fq})
        samtools fastq {input.bam} | gzip > {output.fq} 2> {log}
        """

rule prepare_kallisto_sleuth:
    input:
        old_fastqs=expand(
            "resources/rna_seq_old/fastqs/{accession}.fastq.gz",
            accession=lambda wildcards: config[f"rna_accessions_old"].keys(),
        ),
        new_fastqs=expand(
            "resources/rna_seq_new/fastqs/{accession}.fastq.gz",
            accession=lambda wildcards: config[f"rna_accessions_new"].keys(),
        ),
    output:
        samples="config/samples.tsv",
        units="config/units.tsv",
    conda:
        "../envs/python.yaml"
    params:
        labels_old=lambda wildcards: config[f"rna_accessions_old"],
        labels_new=lambda wildcards: config[f"rna_accessions_new"],
    log:
        "logs/prepare_kallisto_sleuth.log",
    script:
        "../scripts/prepare_kallisto_sleuth.py"

# rule copy_samples_and_units:
#     input:
#         samples="config/samples_{rna_type}.tsv",
#         units="config/units_{rna_type}.tsv",
#     output:
#         samples="{rna_data}/{base}/config/samples_{rna_type}.tsv",
#         units="{rna_data}/{base}/config/units_{rna_type}.tsv",
#     conda:
#         "../envs/python.yaml"
#     log:
#         "logs/copy_samples_units_{rna_data}_{base}_{rna_type}.log",
#     shell:
#         "cp {input.samples} {output.samples} && cp {input.units} {output.units}"


# We need to copy the fastqs from {rna_data}/resources/fastqs to {rna_data}/resources/fastqs
# so that the kallisto-sleuth modules can find them. This problem occurs because
# the kallisto-sleuth modules are prefixed with rna_{rna_data}/base_{base},
# and the fastqs are stored at rna_{rna_data}/resources/fastqs/{sample}.fastq.gz. There is no option to have an output prefix only for modules.
# rule copy_fastqs:
#     input:
#         fq="{rna_data}/resources/fastqs/{sample}.fastq.gz",
#     output:
#         fq="{rna_data}/{base}/resources/fastqs/{sample}.fastq.gz",
#     shell:
#         """
#         cp {input.fq} {output.fq}
#         """




# ---------------------------------------------------------------------------
# Compare differential methylation (DMR) to differential expression
# ---------------------------------------------------------------------------
# The {rna_data} wildcard is always "rna_old" or "rna_new".
# The kallisto-sleuth modules are prefixed with rna_{rna_data}/base_{base},
# so the diffexp tables live at:
#   rna_old/base_{base}/results/tables/diffexp/condition.genes-representative.diffexp_postprocessed.tsv
#   rna_new/base_{base}/results/tables/diffexp/condition.genes-representative.diffexp_postprocessed.tsv


def compute_diffexp_tables(wildcards, directive):
    computed_diffexp = config["kallisto_sleuth"]["diffexp"]["models"].keys()
    results = []
    # Remove the trailing "derm" from every cell_type
    for cell_type in get_non_base_layers(wildcards.base):
        cell_type = cell_type.replace("derm", "")
        rna_type = wildcards.rna_data.replace("rna_", "")
        if f"{cell_type}_vs_{wildcards.base.replace("derm", "")}_{rna_type}" in computed_diffexp:
            results.append(f"results/tables/diffexp/{cell_type}_vs_{wildcards.base.replace("derm", "")}_{rna_type}.genes-representative.diffexp_postprocessed.tsv") if directive == "input" else results.append(1)
        else:
            results.append(f"results/tables/diffexp/{wildcards.base.replace("derm", "")}_vs_{cell_type}_{rna_type}.genes-representative.diffexp_postprocessed.tsv") if directive == "input" else results.append(-1)
    return results



rule compare_dmr_to_diffexp_no_tfs:
    input:
        diffexp = lambda wildcards: compute_diffexp_tables(wildcards, "input"),
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
        val_genes="resources/rna_seq/val_genes.tsv",
    output:
        tsv="results/{platform}/{caller}/base_{base}/{rna_data}/diffexp_vs_dmrs_{annotation_type}.tsv",
        html="results/{platform}/{caller}/base_{base}/{rna_data}/diffexp_vs_dmrs_{annotation_type}.html",
        val_genes="results/{platform}/{caller}/base_{base}/{rna_data}/val_genes_{annotation_type}.tsv",
    # wildcard_constraints:
    #     # rna_data is strictly rna_old or rna_new – no slashes
    #     rna_data="old|new",
    conda:
        "../envs/python.yaml"
    params:
        annotation_type=lambda wildcards: wildcards.annotation_type,
        base=lambda wildcards: wildcards.base,
        diffexp_base_signs = lambda wildcards: compute_diffexp_tables(wildcards, "params"),
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
        val_genes="results/{platform}/{caller}/base_{base}/{rna_data}/val_genes_{annotation_type}.tsv",
    output:
        report(
            directory(
                "results/{platform}/{caller}/base_{base}/{rna_data}/diffexp_vs_dmrs_no_tfs_{annotation_type}"
            ),
            caption="../report/diffexp_vs_dmrs.rst",
            htmlindex="index.html",
            category=lambda wildcards: f"DiffExp-DMRs Comparison - {wildcards.rna_data}",
            subcategory=f"No tfs",

            # subcategory=f"{wildcards.annotation_type}",
            labels=lambda wildcards: {
                "base": wildcards.base,
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
            category=lambda wildcards: f"DiffExp-DMRs Comparison - {wildcards.rna_data}",
            subcategory= f"With tfs",
            labels=lambda wildcards: {
                "base": wildcards.base,
            },
        ),
    wildcard_constraints:
        rna_data="rna_old|rna_new",
    log:
        "logs/diffexp_dmvzrd/diffexp_dmr_datavzrd_with_tfs/{platform}_{caller}_{base}_{rna_data}_{annotation_type}.log",
    wrapper:
        "v9.2.0/utils/datavzrd"


rule view_val_genes:
    input:
        config=workflow.source_path("../resources/val_genes.rds"),
    output:
        tsv="resources/rna_seq/val_genes.tsv"
    conda:
        "../envs/enrichment.yaml"
    log:
        "logs/val_genes/val_genes.log",
    script:
        "../scripts/view_val_genes.R"

# rule val_genes_table:
#     input:
#         "results/val_genes.tsv",
#         "results/{platform}/{caller}/base_{base}/{rna_data}/diffexp_vs_dmrs_{annotation_type}.tsv"
#     output:
#         tsv="results/{platform}/{caller}/base_{base}/{rna_data}/val_genes_{annotation_type}.tsv"
#     conda:
#         "../envs/python.yaml"
#     log:
#         "logs/val_genes/val_genes_{platform}_{caller}_{base}_{rna_data}_{annotation_type}.log",
#     script:
#         "../scripts/filter_val_genes.py"
