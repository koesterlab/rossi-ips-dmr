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

