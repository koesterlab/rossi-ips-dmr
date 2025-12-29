
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
        diffexp="results/tables/diffexp/condition.genes-representative.diffexp_postprocessed.tsv",
        # diffexp="results/tables/diffexp/condition.genes-aggregated.diffexp.tsv",
        endoderm="results/{platform}/{caller}/dmr_calls/endoderm/genes_transcripts/chipseeker_postprocessed_filtered.tsv",
        mesoderm="results/{platform}/{caller}/dmr_calls/mesoderm/genes_transcripts/chipseeker_postprocessed_filtered.tsv",
        ectoderm="results/{platform}/{caller}/dmr_calls/ectoderm/genes_transcripts/chipseeker_postprocessed_filtered.tsv",
    output:
        report(
            "results/{platform}/{caller}/rna_seq_comp/diffexp_vs_dmrs.html",
            caption="../report/rna_seq.rst",
            category="DiffExp-Methylation Comparison",
            subcategory=lambda wildcards: f"{wildcards.platform} - {wildcards.caller}",
        ),
    conda:
        "../envs/python.yaml"
    log:
        "logs/compare_diffexp/{platform}_{caller}.log",
    script:
        "../scripts/compare_diffexp_dmrs.py"
