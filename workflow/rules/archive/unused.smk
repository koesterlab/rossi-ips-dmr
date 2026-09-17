# Rules that are currently not used by the workflow. Kept for reference; not included.

rule download_varlociraptor:
    output:
        directory(
            "resources/tools/varlociraptor",
        ),
    log:
        "logs/download_varlociraptor.log",
    shell:
        """
        PARENT_DIR=$(dirname {output})
        mkdir -p $PARENT_DIR
        cd $PARENT_DIR
        git clone git@github.com:varlociraptor/varlociraptor.git
        cd varlociraptor
        git checkout methylation-paired-end
        """


rule calls_to_vcf:
    input:
        "results/{platform}/varlo/meth_calling/{sample}/{fdr}/calls_{scatteritem}.filtered.bcf",
    output:
        "results/{platform}/varlo/meth_calling/{sample}/{fdr}/calls_{scatteritem}.vcf",
    conda:
        "../envs/samtools.yaml"
    threads: 10
    log:
        "logs/varlociraptor/calls_to_vcf/{platform}_{sample}_{scatteritem}_{fdr}.log"
    shell:
        """
        bcftools view --threads {threads} {input} -o {output} 2> {log}
        """

rule mosdepth:
    input:
        bam="resources/{platform}/{sample}.bam",
        bai="resources/{platform}/{sample}.bam.bai",
        bed="resources/candidates/candidates.bed"
    output:
        dist="resources/{platform}/{sample}.mosdepth.global.dist.txt",
        per_base="resources/{platform}/{sample}.per-base.bed.gz",
        summary="resources/{platform}/{sample}.mosdepth.summary.txt",
    log:
        "logs/mosdepth/{platform}_{sample}.log",
    params:
        extra="--fast-mode",
    threads: 4
    wrapper:
        "v9.4.2/bio/mosdepth"



rule merge_mosdepth:
    input:
        expand(
            "resources/{{platform}}/{sample}.per-base.bed.gz",
            sample=ALL_GERM_LAYERS
        )
    output:
        "resources/{platform}/coverage.parquet",
    log:
        "logs/merge_mosdepth/{platform}.log",
    script:
        "../../scripts/archive/merge_mosdepth.py"

rule scatter_plot:
    input:
        calls="results/{platform}/{caller}/meth_calling/calls_0.05.parquet",
    output:
        report(
            "results/{platform}/{caller}/base_{base}/plots_paper/{group2}/scatter_plot.png",
            caption="../report/scatter_plot.rst",
            category="Plots paper",
            subcategory=lambda wildcards: f"{wildcards.platform} - {wildcards.caller}",
            labels=lambda wildcards: {
                "Plot": "1B",
                "Base": wildcards.base,
                "Type": wildcards.group2,
            },
        ),
    params:
        group1=lambda wildcards: wildcards.base,
        group2=lambda wildcards: wildcards.group2,
    resources:
        mem_mb=16000,
    conda:
        "../envs/python.yaml"
    log:
        "logs/scatter_plot/{platform}_{caller}_{base}_{group2}.log",
    script:
        "../scripts/scatter_plot.py"


rule scatter_plot_endo_meso:
    input:
        calls="results/{platform}/{caller}/meth_calling/calls_1.0.parquet",
    output:
        report(
            "results/{platform}/{caller}/plots_paper/endo_meso/scatter_plot.png",
            caption="../report/scatter_plot.rst",
            category="Plots paper",
            subcategory=lambda wildcards: f"{wildcards.platform} - {wildcards.caller}",
            labels=lambda wildcards: {
                "Plot": "1C",
                "Type": "endo_meso",
            },
        ),
    resources:
        mem_mb=16000,
    params:
        group1="mesoderm",
        group2="endoderm",
    conda:
        "../envs/python.yaml"
    log:
        "logs/scatter_plot_endo_meso/{platform}_{caller}.log",
    script:
        "../scripts/scatter_plot.py"

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
