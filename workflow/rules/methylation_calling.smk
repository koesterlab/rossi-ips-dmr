
scattergather:
    split_candidates=config["scatter_items"],


rule find_candidates:
    input:
        fasta="resources/genome.fasta",
    output:
        "resources/candidates/candidates.bcf",
    log:
        "logs/find_candidates.log",
    conda:
        "../envs/varlociraptor.yaml"
    shell:
        """ 
        varlociraptor -- methylation-candidates {input.fasta} {output} --motifs CG
        """


rule split_candidates:
    input:
        "resources/candidates/candidates.bcf",
    output:
        scatter.split_candidates("resources/candidates/candidates_{scatteritem}.bcf"),
    log:
        "logs/split_candidates.log",
    conda:
        "../envs/rbt.yaml"
    shell:
        "rbt vcf-split {input} {output}"


# We want to extract only the important reads from the bam file. Since the candidates span more than one chromosome we cant just use the normal samtools view -b
# rule candidates_to_bed:
#     input:
#         "resources/candidates/candidates_{scatteritem}.bcf",
#     output:
#         "resources/candidates/candidates_{scatteritem}.bed",
#     log:
#         "logs/candidates_to_bed_{scatteritem}.log",
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
    params:
        sequencer=lambda wildcards: samples[wildcards.sample][1],
    shell:
        """ 
        varlociraptor -- preprocess variants {input.genome} --candidates {input.candidates} --bam {input.alignments} --read-type {params.sequencer} --max-depth 1000 > {output}
        """


# rule call_methylation_single:
# input:
#     varlo_path="resources/tools/varlociraptor",
#     preprocess_obs="results/{platform}/varlo/meth_calling/{sample}/normal_{scatteritem}.bcf",
#     scenario="resources/scenario.yaml",
# output:
#     "results/{platform}/varlo/meth_calling/{sample}/calls_{scatteritem}.bcf",
# conda:
#     "../envs/varlociraptor.yaml"
# shell:
#     """
#     cd {input.varlo_path}
#     cargo run --release -- call variants --omit-strand-bias generic --scenario {input.scenario} --obs normal={input.preprocess_obs} > {output}
#     """


# rule call_methylation_together:
#     input:
#         varlo_path="resources/tools/varlociraptor",
#         pb_psc="results/pacbio/varlo/meth_calling/psc/normal_{scatteritem}.bcf",
#         pb_ecto="results/pacbio/varlo/meth_calling/ectoderm/normal_{scatteritem}.bcf",
#         pb_endo="results/pacbio/varlo/meth_calling/endoderm/normal_{scatteritem}.bcf",
#         pb_meso="results/pacbio/varlo/meth_calling/mesoderm/normal_{scatteritem}.bcf",
#         np_psc="results/nanopore/varlo/meth_calling/psc/normal_{scatteritem}.bcf",
#         np_ecto="results/nanopore/varlo/meth_calling/ectoderm/normal_{scatteritem}.bcf",
#         np_endo="results/nanopore/varlo/meth_calling/endoderm/normal_{scatteritem}.bcf",
#         np_meso="results/nanopore/varlo/meth_calling/mesoderm/normal_{scatteritem}.bcf",
#         scenario="resources/scenario.yaml",
#     output:
#         "results/platforms_combined/varlo/meth_calling/cell_lines_combined/calls_{scatteritem}.bcf",
#     conda:
#         "../envs/varlociraptor.yaml"
#     log:
#         "logs/compute_meth_together_{scatteritem}.log",
#     resources:
#         mem_mb=128000,
#     shell:
#         """ 
#         cd {input.varlo_path}
#         cargo run --release -- call variants --omit-strand-bias generic --scenario {input.scenario} \
#             --obs psc_pacbio={input.pb_psc} meso_pacbio={input.pb_meso} endo_pacbio={input.pb_endo} ecto_pacbio={input.pb_ecto} \
#             psc_nanopore={input.np_psc} meso_nanopore={input.np_meso} endo_nanopore={input.np_endo} ecto_nanopore={input.np_ecto}  > {output} 2> {log}
#         """

rule call_methylation:
    input:
        pb="results/pacbio/varlo/meth_calling/{sample}/normal_{scatteritem}.bcf",
        np="results/nanopore/varlo/meth_calling/{sample}/normal_{scatteritem}.bcf",
        scenario="resources/scenario_common.yaml",
    output:
        "results/platforms_combined/varlo/meth_calling/{sample}/calls_{scatteritem}.bcf",
    conda:
        "../envs/varlociraptor.yaml"
    log:
        "logs/compute_meth_together_{sample}_{scatteritem}.log",
    resources:
        mem_mb=128000,
    shell:
        """ 
        varlociraptor -- call variants generic --scenario {input.scenario} \
            --obs pacbio={input.pb}  nanopore={input.np}   > {output} 2> {log}
        """



# # TODO: Reactivate, right now it deletes too much data
# rule filter_calls:
#     input:
#         calls="results/{platform}/varlo/meth_calling/{sample}/calls_{scatteritem}.bcf",
#         varlo_path="resources/tools/varlociraptor",
#     output:
#         "results/{platform}/varlo/meth_calling/{sample}/calls_{scatteritem}.filtered.bcf",
#     conda:
#         "../envs/varlociraptor.yaml"
#     params:
#         event="PRESENT",
#     shell:
#         """
#         PIPELINE_PATH=$(pwd)
#         cd {input.varlo_path}
#         cargo run --release -- filter-calls control-fdr --mode local-smart {input.calls} --events {params.event} --fdr 1 > {output}
#         """


rule calls_to_vcf:
    input:
        "results/{platform}/varlo/meth_calling/{sample}/calls_{scatteritem}.bcf",
    output:
        "results/{platform}/varlo/meth_calling/{sample}/calls_{scatteritem}.vcf",
    conda:
        "../envs/samtools.yaml"
    threads: 10
    shell:
        """
        bcftools view --threads {threads} {input} -o {output}
        """


rule gather_calls:
    input:
        gather.split_candidates(
            "results/{{platform}}/varlo/meth_calling/{{sample}}/calls_{scatteritem}.vcf"
        ),
    output:
        "results/{platform}/varlo/meth_calling/{sample}/varlo.vcf",
    conda:
        "../envs/cat.yaml"
    shell:
        "cat {input} > {output}"


rule df_from_calls:
    input:
        undifferentiated="results/{platform}/{caller}/meth_calling/psc/{caller}.vcf",
        meso="results/{platform}/{caller}/meth_calling/mesoderm/{caller}.vcf",
        endo="results/{platform}/{caller}/meth_calling/endoderm/{caller}.vcf",
        ecto="results/{platform}/{caller}/meth_calling/ectoderm/{caller}.vcf",
        # "results/{platform}/{caller}/meth_calling/{group}/calls.vcf",
    output:
        "results/{platform}/{caller}/meth_calling/calls.parquet",
    conda:
        "../envs/plot.yaml"
    log:
        "logs/df_from_calls/{platform}/{caller}.log"
    resources:
        mem_mb=64000,
    params:
        meth_caller=lambda wildcards: wildcards.caller,
        prob_pres_threshhold=config["prob_pres_threshold"],
        prob_abs_threshhold=config["prob_abs_threshold"],
        alpha=config["alpha"],
    script:
        "../scripts/df_from_calls.py"