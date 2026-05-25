
scattergather:
    split_candidates=config["scatter_items"],


rule find_candidates:
    input:
        fasta="resources/genome.fasta",
    output:
        "resources/candidates/candidates.bcf",
    log:
        "logs/varlociraptor/find_candidates.log",
    conda:
        "../envs/varlociraptor.yaml"
    shell:
        """
        varlociraptor methylation-candidates {input.fasta} {output} --motifs CG 2> {log}
        """


rule split_candidates:
    input:
        "resources/candidates/candidates.bcf",
    output:
        scatter.split_candidates("resources/candidates/candidates_{scatteritem}.bcf"),
    log:
        "logs/varlociraptor/split_candidates.log",
    conda:
        "../envs/rbt.yaml"
    shell:
        "rbt vcf-split {input} {output} 2> {log}"


# # We want to extract only the important reads from the bam file. Since the candidates span more than one chromosome we cant just use the normal samtools view -b
# rule candidates_to_bed:
#     input:
#         "resources/candidates/candidates_{scatteritem}.bcf",
    # output:
#         "resources/candidates/candidates_{scatteritem}.bed",
#     log:
#         "logs/varlociraptor/candidates_to_bed/{scatteritem}.log",
#     conda:
#         "../envs/pysam.yaml"
#     script:
#         "../scripts/candidates_to_bed.py"


rule compute_meth_observations:
    input:
        genome="resources/genome.fasta",
        genomeIndex="resources/genome.fasta.fai",
        # alignments="resources/{platform}/meth_calling/{sample}.bam",
        alignments="resources/{platform}/{sample}/alignment_{scatteritem}.bam",
        alignments_index="resources/{platform}/{sample}/alignment_{scatteritem}.bam.bai",
        # alignments_index="resources/{platform}/meth_calling/{sample}.bam.bai",
        candidates="resources/candidates/candidates_{scatteritem}.bcf",
    output:
        "results/{platform}/varlo/meth_calling/{sample}/normal_{scatteritem}.bcf",
    conda:
        "../envs/varlociraptor.yaml"
    log:
        "logs/varlociraptor/compute_meth_observations/{platform}_{sample}_{scatteritem}.log",
    shell:
        """
        varlociraptor preprocess variants {input.genome} --candidates {input.candidates} --bam {input.alignments}  --methylation-read-type annotated --max-depth 1000 > {output} 2> {log}
        """


rule call_methylation_single:
    input:
        preprocess_obs="results/{platform}/varlo/meth_calling/{sample}/normal_{scatteritem}.bcf",
        scenario=workflow.source_path("../resources/scenarios/scenario.yaml"),
    output:
        "results/{platform}/varlo/meth_calling/{sample}/calls_{scatteritem}.bcf",
    conda:
        "../envs/varlociraptor.yaml"
    wildcard_constraints:
        platform="(pacbio|nanopore)",
    log:
        "results/call_methylation_single/{platform}_{sample}_{scatteritem}.log",
    shell:
        """
            varlociraptor call variants generic --scenario {input.scenario} --obs normal={input.preprocess_obs} > {output} 2> {log}
        """

rule call_methylation:
    input:
        pb="results/pacbio/varlo/meth_calling/{sample}/normal_{scatteritem}.bcf",
        np="results/nanopore/varlo/meth_calling/{sample}/normal_{scatteritem}.bcf",
        scenario=workflow.source_path("../resources/scenarios/scenario_common.yaml"),
    output:
        "results/platforms_combined/varlo/meth_calling/{sample}/calls_{scatteritem}.bcf",
    conda:
        "../envs/varlociraptor.yaml"
    log:
        "logs/varlociraptor/compute_meth_together/{sample}_{scatteritem}.log",
    resources:
        mem_mb=128000,
    shell:
        """
        varlociraptor call variants generic --scenario {input.scenario} \
            --obs pacbio={input.pb}  nanopore={input.np}   > {output} 2> {log}
        """

rule filter_prob_absent:
    input:
        "results/{platform}/varlo/meth_calling/{sample}/calls_{scatteritem}.bcf",
    output:
        "results/{platform}/varlo/meth_calling/{sample}/calls_{scatteritem}.filtered_absent.bcf",
    conda:
        "../envs/varlociraptor.yaml"
    params:
        event="ABSENT",
        fdr=config["fdr_absent"]
    log:
        "logs/varlociraptor/filter_prob_absent/{platform}_{sample}_{scatteritem}.log",
    shell:
        "varlociraptor filter-calls control-fdr --mode local-smart {input} --events {params.event} --fdr {params.fdr} > {output} 2> {log}"

rule filter_prob_present:
    input:
        "results/{platform}/varlo/meth_calling/{sample}/calls_{scatteritem}.bcf",
    output:
        "results/{platform}/varlo/meth_calling/{sample}/calls_{scatteritem}.filtered_present.bcf",
    conda:
        "../envs/varlociraptor.yaml"
    params:
        event="PRESENT",
        fdr=config["fdr_present"]
    log:
        "logs/varlociraptor/filter_prob_present/{platform}_{sample}_{scatteritem}.log",
    shell:
        "varlociraptor filter-calls control-fdr --mode local-smart {input} --events {params.event} --fdr {params.fdr} > {output} 2> {log}"

rule concatenate_filtered_calls:
    input:
        absent="results/{platform}/varlo/meth_calling/{sample}/calls_{scatteritem}.filtered_absent.bcf",
        absent_index="results/{platform}/varlo/meth_calling/{sample}/calls_{scatteritem}.filtered_absent.bcf.csi",
        present="results/{platform}/varlo/meth_calling/{sample}/calls_{scatteritem}.filtered_present.bcf",
        present_index="results/{platform}/varlo/meth_calling/{sample}/calls_{scatteritem}.filtered_present.bcf.csi",
    output:
        "results/{platform}/varlo/meth_calling/{sample}/calls_{scatteritem}.filtered.bcf",
    log:
        "logs/varlociraptor/concatenate_filtered_calls/{platform}_{sample}_{scatteritem}.log",
    shell:
        "bcftools concat -a {input.absent} {input.present} -o {output} 2> {log}"



rule calls_to_vcf:
    input:
        "results/{platform}/varlo/meth_calling/{sample}/calls_{scatteritem}.filtered.bcf",
    output:
        "results/{platform}/varlo/meth_calling/{sample}/calls_{scatteritem}.vcf",
    conda:
        "../envs/samtools.yaml"
    threads: 10
    log:
        "logs/varlociraptor/calls_to_vcf/{platform}_{sample}_{scatteritem}.log"
    shell:
        """
        bcftools view --threads {threads} {input} -o {output} 2> {log}
        """


rule gather_calls:
    input:
        gather.split_candidates(
            "results/{{platform}}/varlo/meth_calling/{{sample}}/calls_{scatteritem}.filtered.bcf"
        ),
    output:
        "results/{platform}/varlo/meth_calling/{sample}/varlo.bcf",
    conda:
        "../envs/samtools.yaml"
    log:
        "logs/varlociraptor/gather_calls/{platform}_{sample}.log",
    shell:
        """
        bcftools concat  {input} -o {output} 2> {log}
        """

rule prepare_wsabi:
    input:
        "results/{platform}/varlo/meth_calling/{sample}/varlo.bcf",
    output:
        "results/wsabi/{platform}_{sample}.tsv.gz",
    conda:
        "../envs/pysam.yaml"
    log:
        "logs/varlociraptor/prepare_wasabi/{platform}_{sample}.log"
    script:
        "../scripts/bcf_to_wasabi_tsv.py"


rule index_bcf:
    input:
        "{bcf}.bcf"
    output:
        "{bcf}.bcf.csi"
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        bcftools index -c {input}
        """


rule df_from_calls:
    input:
        undifferentiated="results/{platform}/{caller}/meth_calling/psc/{caller}.bcf",
        undifferentiated_index="results/{platform}/{caller}/meth_calling/psc/{caller}.bcf.csi",
        meso="results/{platform}/{caller}/meth_calling/mesoderm/{caller}.bcf",
        meso_index="results/{platform}/{caller}/meth_calling/mesoderm/{caller}.bcf.csi",
        endo="results/{platform}/{caller}/meth_calling/endoderm/{caller}.bcf",
        endo_index="results/{platform}/{caller}/meth_calling/endoderm/{caller}.bcf.csi",
        ecto="results/{platform}/{caller}/meth_calling/ectoderm/{caller}.bcf",
        ecto_index="results/{platform}/{caller}/meth_calling/ectoderm/{caller}.bcf.csi",
        # "results/{platform}/{caller}/meth_calling/{group}/calls.vcf",
    output:
        "results/{platform}/{caller}/meth_calling/calls.parquet",
    conda:
        "../envs/pysam.yaml"
    log:
        "logs/varlociraptor/df_from_calls/{platform}_{caller}.log"
    resources:
        mem_mb=64000,
    params:
        meth_caller=lambda wildcards: wildcards.caller,
    script:
        "../scripts/df_from_calls.py"
