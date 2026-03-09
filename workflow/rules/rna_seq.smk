
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
        samples="rna_old/config/old_samples.tsv",
        units="rna_old/config/old_units.tsv",
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
        expand("rna_new/resources/fastqs/{barcode}", barcode=[f"SQK-NBD114-24_barcode{str(i).zfill(2)}.bam" for i in range(1, 13)]),
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
        fq="rna_new/resources/fastqs/{sample}.fastq.gz"
    conda:
        "../envs/samtools.yaml"
    log:
        "logs/bam_to_fastq/{sample}.log"
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
        samples="rna_new/config/new_samples.tsv",
        units="rna_new/config/new_units.tsv",
    conda:
        "../envs/python.yaml"
    log:
        "logs/prepare_kallisto_sleuth_new.log",
    params:
        types=config["rna_accessions_new"],
    script:
        "../scripts/prepare_kallisto_sleuth_new.py"



# Compare the differential expression results with the DMR associated genes
rule compare_dmr_to_diffexp_no_tfs:
    input:
        diffexp="{rna_data}/results/tables/diffexp/condition.genes-representative.diffexp_postprocessed.tsv",
        # diffexp="results/tables/diffexp/condition.genes-aggregated.diffexp.tsv",
        endoderm="results/{platform}/{caller}/dmr_calls/endoderm/genes_transcripts/chipseeker_postprocessed.tsv",
        mesoderm="results/{platform}/{caller}/dmr_calls/mesoderm/genes_transcripts/chipseeker_postprocessed.tsv",
        ectoderm="results/{platform}/{caller}/dmr_calls/ectoderm/genes_transcripts/chipseeker_postprocessed.tsv",
    output:
        "results/{platform}/{caller}/{rna_data}/diffexp_vs_dmrs_{annotation_type}.tsv",
        "results/{platform}/{caller}/{rna_data}/diffexp_vs_dmrs_{annotation_type}.html",
    conda:
        "../envs/python.yaml"
    params:
        annotation_type=lambda wildcards: wildcards.annotation_type,
    log:
        "logs/compare_dmr_to_diffexp_no_tfs/{platform}_{caller}_{annotation_type}_{rna_data}.log",
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
        comp="results/{platform}/{caller}/{rna_data}/diffexp_vs_dmrs_{annotation_type}.tsv",
    output:
        focus_tfs="results/{platform}/{caller}/{rna_data}/tfs/diffexp_vs_dmrs_tfs_only_{annotation_type}.tsv",
        comp_tf_adj="results/{platform}/{caller}/{rna_data}/tfs/diffexp_vs_dmrs_{annotation_type}.tsv",
        plot="results/{platform}/{caller}/{rna_data}/tfs/diffexp_vs_dmrs_{annotation_type}.html",
    conda:
        "../envs/python.yaml"
    params:
        annotation_type=lambda wildcards: wildcards.annotation_type,
    log:
        "logs/compare_dmr_to_diffexp_with_tfs/{platform}_{caller}_{annotation_type}_{rna_data}_tfs.log",
    script:
        "../scripts/compare_dmr_to_diffexp_with_tfs.py"


# We can't use the datavzrd wrapper since we need the latest datvzrd version.
rule datavzrd_dmr_vs_diffexp_no_tfs:
    input:
        config=workflow.source_path("../resources/dmr_vs_diffexp_no_tfs.yaml"),
        table="results/{platform}/{caller}/{rna_data}/diffexp_vs_dmrs_{annotation_type}.tsv",
    output:
        report(
            directory(
                "results/{platform}/{caller}/{rna_data}/diffexp_vs_dmrs_no_tfs_{annotation_type}"
            ),
            caption="../report/diffexp_vs_dmrs.rst",
            htmlindex="index.html",
            category="DiffExp-DMRs Comparison",
            subcategory=lambda wildcards: f"Comparisons",

            # subcategory=lambda wildcards: f"{wildcards.annotation_type}",
            labels=lambda wildcards: {
                "type": "no transcription factors",
            },
        ),
    log:
        "logs/diffexp_dmvzrd/diffexp_dmr_datavzrd/{platform}_{caller}_{annotation_type}_{rna_data}.log",
    conda:
        "../envs/datavzrd.yaml"
    wrapper:
        # "641c90c4da86d4acf2022f347f3c8017334c0f44/utils/datavzrd"
        "v9.2.0/utils/datavzrd"


rule datavzrd_dmr_vs_diffexp_with_tfs:
    input:
        config=workflow.source_path("../resources/dmr_vs_diffexp_with_tfs.yaml"),
        complete="results/{platform}/{caller}/{rna_data}/tfs/diffexp_vs_dmrs_{annotation_type}.tsv",
        focus_tfs="results/{platform}/{caller}/{rna_data}/tfs/diffexp_vs_dmrs_tfs_only_{annotation_type}.tsv",
    output:
        report(
            directory(
                "results/{platform}/{caller}/{rna_data}/diffexp_vs_dmrs_with_tfs_{annotation_type}"
            ),
            caption="../report/diffexp_vs_dmrs.rst",
            htmlindex="index.html",
            category="DiffExp-DMRs Comparison",
            subcategory=lambda wildcards: f"Comparisons",
            labels=lambda wildcards: {
                "type": "with transcription factors",
            },
        ),
    log:
        "logs/diffexp_dmvzrd/diffexp_dmr_datavzrd/{platform}_{caller}_{annotation_type}_{rna_data}.log",
    wrapper:
        # "641c90c4da86d4acf2022f347f3c8017334c0f44/utils/datavzrd"
        "v9.2.0/utils/datavzrd"
