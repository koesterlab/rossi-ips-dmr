
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
rule compare_diffexp:
    input:
        diffexp="results/tables/diffexp/condition.genes-representative.diffexp.tsv",
        # diffexp="results/tables/diffexp/condition.genes-aggregated.diffexp.tsv",
        endoderm="results/{platform}/{caller}/dmr_calls/endoderm/genes_transcripts/chipseeker_postprocessed_complete.tsv",
        mesoderm="results/{platform}/{caller}/dmr_calls/mesoderm/genes_transcripts/chipseeker_postprocessed_complete.tsv",
        ectoderm="results/{platform}/{caller}/dmr_calls/ectoderm/genes_transcripts/chipseeker_postprocessed_complete.tsv",
    output:
        report(
            "results/{platform}/{caller}/rna_seq_comp/diffexp_vs_dmrs_{annotation_type}.html",
            caption="../report/rna_seq.rst",
            category="DiffExp-Methylation Comparison",
            subcategory=lambda wildcards: f"{wildcards.annotation_type} - scatter",
        ),
        "results/{platform}/{caller}/rna_seq_comp/diffexp_vs_dmrs_{annotation_type}.tsv",
    conda:
        "../envs/python.yaml"
    params:
        annotation_type=lambda wildcards: wildcards.annotation_type,
    log:
        "logs/compare_diffexp/{platform}_{caller}_{annotation_type}.log",
    script:
        "../scripts/compare_diffexp_dmrs.py"


rule diffexp_dmr_datavzrd:
    input:
        config=workflow.source_path("../resources/diffexp_vs_dmrs.yaml"),
        comp="results/{platform}/{caller}/rna_seq_comp/tfs/diffexp_vs_dmrs_promoter.tsv",
    output:
        report(
            directory(
                "results/{platform}/{caller}/rna_seq_comp/diffexp_vs_dmrs_promoter"
            ),
            caption="../report/diffexp_vs_dmrs.rst",
            htmlindex="index.html",
            category="DiffExp-Methylation Comparison",
            subcategory=lambda wildcards: f"{wildcards.annotation_type} - table",
        ),
    log:
        "logs/diffexp_dmvzrd/diffexp_dmr_datavzrd/{platform}_{caller}.log",
    wrapper:
        "v3.13.0/utils/datavzrd"


# Compare the differential expression results with the DMR associated genes
rule compute_tfs:
    input:
        diffexp="results/tables/diffexp/condition.genes-representative.diffexp.tsv",
        samples="config/samples.tsv",
    output:
        # report(
        #     "results/{platform}/{caller}/rna_seq_comp/tfs/diffexp_vs_dmrs_{annotation_type}.html",
        #     caption="../report/rna_seq.rst",
        #     category="DiffExp-Methylation Comparison",
        #     subcategory=lambda wildcards: f"{wildcards.annotation_type}",
        # ),
        "results/{platform}/{caller}/rna_seq_comp/tfs/{layer}_all_effects.tsv",
        "results/{platform}/{caller}/rna_seq_comp/tfs/{layer}_significant_tfs.tsv",
    conda:
        "../envs/decoupler.yaml"
    params:
        layer=lambda wildcards: wildcards.layer,
        max_p_value=0.05,
        minsize=5,
        min_contribution=1,
    log:
        "logs/compare_diffexp/{platform}_{caller}_{layer}.log",
    script:
        "../scripts/compute_tfs_target_genes.R"


# Compare the differential expression results with the DMR associated genes
rule compare_diffexp_with_tfs:
    input:
        diffexp="results/tables/diffexp/condition.genes-representative.diffexp.tsv",
        # diffexp="results/tables/diffexp/condition.genes-aggregated.diffexp.tsv",
        endoderm="results/{platform}/{caller}/dmr_calls/endoderm/genes_transcripts/chipseeker_postprocessed_complete.tsv",
        mesoderm="results/{platform}/{caller}/dmr_calls/mesoderm/genes_transcripts/chipseeker_postprocessed_complete.tsv",
        ectoderm="results/{platform}/{caller}/dmr_calls/ectoderm/genes_transcripts/chipseeker_postprocessed_complete.tsv",
        ectoderm_tfs="results/{platform}/{caller}/rna_seq_comp/tfs/ectoderm_significant_tfs.tsv",
        mesoderm_tfs="results/{platform}/{caller}/rna_seq_comp/tfs/mesoderm_significant_tfs.tsv",
        endoderm_tfs="results/{platform}/{caller}/rna_seq_comp/tfs/endoderm_significant_tfs.tsv",
    output:
        report(
            "results/{platform}/{caller}/rna_seq_comp/tfs/diffexp_vs_dmrs_{annotation_type}.html",
            caption="../report/rna_seq.rst",
            category="DiffExp-Methylation Comparison",
            subcategory=lambda wildcards: f"{wildcards.annotation_type} - scatter tfs",
        ),
        "results/{platform}/{caller}/rna_seq_comp/tfs/diffexp_vs_dmrs_{annotation_type}.tsv",
    conda:
        "../envs/python.yaml"
    params:
        annotation_type=lambda wildcards: wildcards.annotation_type,
    log:
        "logs/compare_diffexp/{platform}_{caller}_{annotation_type}.log",
    script:
        "../scripts/compare_diffexp_dmrs_with_tfs.py"
