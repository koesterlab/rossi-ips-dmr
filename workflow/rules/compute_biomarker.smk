
rule candidates_to_bed:
    input:
        "resources/candidates/{candidates}.bcf",
    output:
        "resources/candidates/{candidates}.bed",
    conda:
        "../envs/samtools.yaml"
    log:
        "logs/mason/candidates_to_bed/{candidates}.log",
    shell:
        """
        bcftools query -f '%CHROM\t%POS\t%REF\n' {input} 2> {log} | \
        awk '{{print $1 "\t" $2-1 "\t" $2-1+length($3)}}' > {output}
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
        "../scripts/merge_mosdepth.py"


rule compute_own_biomarker:
    input:
        nanopore="results/nanopore/{caller}/meth_calling/calls_0.01.parquet",
        pacbio="results/pacbio/{caller}/meth_calling/calls_0.01.parquet",
        coverage_nanopore="resources/nanopore/coverage.parquet",
        coverage_pacbio="resources/pacbio/coverage.parquet",
    output:
        "results/platforms_combined/{caller}/plots_paper/own_biomarker.parquet"
    resources:
        mem_mb=16000,
    conda:
        "../envs/python.yaml"
    log:
        "logs/compute_own_biomarker/{caller}.log",
    script:
        "../scripts/compute_own_biomarker.py"
