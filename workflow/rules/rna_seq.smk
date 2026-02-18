
rule get_fastq_se_gz:
    output:
        "resources/rna_seq/fastqs/{accession}.fastq.gz",
    log:
        "logs/get_fastq_se_gz/{accession}.gz.log",
    params:
        extra="--skip-technical",
    threads: 6
    wrapper:
        "v7.6.0/bio/sra-tools/fasterq-dump"


rule prepare_kallisto_sleuth:
    input:
        fastqs=expand(
            "resources/rna_seq/fastqs/{accession}.fastq.gz",
            accession=config["rna_accessions"].keys(),
        ),
    output:
        samples="config/samples.tsv",
        units="config/units.tsv",
    conda:
        "../envs/python.yaml"
    log:
        "logs/prepare_kallisto_sleuth.log",
    params:
        types=config["rna_accessions"],
    script:
        "../scripts/prepare_kallisto_sleuth.py"


# Compare the differential expression results with the DMR associated genes
rule compare_dmr_to_diffexp_no_tfs:
    input:
        diffexp="results/tables/diffexp/condition.genes-representative.diffexp.tsv",
        # diffexp="results/tables/diffexp/condition.genes-aggregated.diffexp.tsv",
        endoderm="results/{platform}/{caller}/dmr_calls/endoderm/genes_transcripts/chipseeker_postprocessed_complete.tsv",
        mesoderm="results/{platform}/{caller}/dmr_calls/mesoderm/genes_transcripts/chipseeker_postprocessed_complete.tsv",
        ectoderm="results/{platform}/{caller}/dmr_calls/ectoderm/genes_transcripts/chipseeker_postprocessed_complete.tsv",
    output:
        "results/{platform}/{caller}/rna_seq_comp/diffexp_vs_dmrs_{annotation_type}.tsv",
        report(
            "results/{platform}/{caller}/rna_seq_comp/diffexp_vs_dmrs_{annotation_type}.html",
            caption="../report/rna_seq.rst",
            category="DiffExp-Methylation Comparison",
            subcategory="scatter plots",
            labels=lambda wildcards: {
                "type": "no transcription factors",
                "annotation type": wildcards.annotation_type,
            },
        ),
    conda:
        "../envs/python.yaml"
    params:
        annotation_type=lambda wildcards: wildcards.annotation_type,
    log:
        "logs/compare_dmr_to_diffexp_no_tfs/{platform}_{caller}_{annotation_type}.log",
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
        comp="results/{platform}/{caller}/rna_seq_comp/diffexp_vs_dmrs_{annotation_type}.tsv",
    output:
        focus_tfs="results/{platform}/{caller}/rna_seq_comp/tfs/diffexp_vs_dmrs_tfs_only_{annotation_type}.tsv",
        comp_tf_adj="results/{platform}/{caller}/rna_seq_comp/tfs/diffexp_vs_dmrs_{annotation_type}.tsv",
        plot="results/{platform}/{caller}/rna_seq_comp/tfs/diffexp_vs_dmrs_{annotation_type}.html",
    conda:
        "../envs/python.yaml"
    params:
        annotation_type=lambda wildcards: wildcards.annotation_type,
    log:
        "logs/compare_dmr_to_diffexp_with_tfs/{platform}_{caller}_{annotation_type}_tfs.log",
    script:
        "../scripts/compare_dmr_to_diffexp_with_tfs.py"


rule datavzrd_dmr_vs_diffexp_no_tfs:
    input:
        config=workflow.source_path("../resources/dmr_vs_diffexp_no_tfs.yaml"),
        table="results/{platform}/{caller}/rna_seq_comp/diffexp_vs_dmrs_{annotation_type}.tsv",
    output:
        report(
            directory(
                "results/{platform}/{caller}/rna_seq_comp/diffexp_vs_dmrs_no_tfs_{annotation_type}"
            ),
            caption="../report/diffexp_vs_dmrs.rst",
            htmlindex="index.html",
            category="DiffExp-DMRs Comparison",
            subcategory=lambda wildcards: f"{wildcards.annotation_type}",
            labels=lambda wildcards: {
                "type": "no transcription factors",
            },
        ),
    log:
        "logs/diffexp_dmvzrd/diffexp_dmr_datavzrd/{platform}_{caller}_{annotation_type}.log",
    wrapper:
        "v9.0.0/utils/datavzrd"


rule datavzrd_dmr_vs_diffexp_with_tfs:
    input:
        config=workflow.source_path("../resources/dmr_vs_diffexp_with_tfs.yaml"),
        complete="results/{platform}/{caller}/rna_seq_comp/tfs/tf_adjusted_{annotation_type}.tsv",
        focus_tfs="results/{platform}/{caller}/rna_seq_comp/tfs/focus_tfs_{annotation_type}.tsv",
    output:
        report(
            directory(
                "results/{platform}/{caller}/rna_seq_comp/diffexp_vs_dmrs_with_tfs_{annotation_type}"
            ),
            caption="../report/diffexp_vs_dmrs.rst",
            htmlindex="index.html",
            category="DiffExp-DMRs Comparison",
            subcategory=lambda wildcards: f"{wildcards.annotation_type}",
            labels=lambda wildcards: {
                "type": "with transcription factors",
            },
        ),
    log:
        "logs/diffexp_dmvzrd/diffexp_dmr_datavzrd/{platform}_{caller}_{annotation_type}.log",
    wrapper:
        "v9.0.0/utils/datavzrd"
