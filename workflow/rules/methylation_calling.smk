
scattergather:
    split_candidates=config["scatter_items"],


rule find_candidates:
    input:
        fasta="resources/genome.fasta",
        varlo_path="resources/tools/varlociraptor",
    output:
        "resources/candidates/candidates.bcf",
    log:
        "logs/find_candidates.log",
    conda:
        "../envs/varlociraptor.yaml"
    shell:
        """ 
        cd {input.varlo_path}
        cargo run -- methylation-candidates {input.fasta} {output}
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
rule candidates_to_bed:
    input:
        "resources/candidates/candidates_{scatteritem}.bcf",
    output:
        "resources/candidates/candidates_{scatteritem}.bed",
    log:
        "logs/candidates_to_bed_{scatteritem}.log",
    conda:
        "../envs/pysam.yaml"
    script:
        "../scripts/candidates_to_bed.py"


rule compute_meth_observations:
    input:
        varlo_path="resources/tools/varlociraptor",
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
        cd {input.varlo_path}
        cargo run --release -- preprocess variants {input.genome} --candidates {input.candidates} --bam {input.alignments} --read-type {params.sequencer} --max-depth 1000 > {output}
        """


rule call_methylation:
    input:
        varlo_path="resources/tools/varlociraptor",
        preprocess_obs="results/{platform}/varlo/meth_calling/{sample}/normal_{scatteritem}.bcf",
        scenario="resources/scenario.yaml",
    output:
        "results/{platform}/varlo/meth_calling/{sample}/calls_{scatteritem}.bcf",
    conda:
        "../envs/varlociraptor.yaml"
    shell:
        """ 
        cd {input.varlo_path}
        cargo run --release -- call variants --omit-strand-bias generic --scenario {input.scenario} --obs normal={input.preprocess_obs} > {output}
        """


# TODO: Reactivate, right now it deletes too much data
rule filter_calls:
    input:
        calls="results/{platform}/varlo/meth_calling/{sample}/calls_{scatteritem}.bcf",
        varlo_path="resources/tools/varlociraptor",
    output:
        "results/{platform}/varlo/meth_calling/{sample}/calls_{scatteritem}.filtered.bcf",
    conda:
        "../envs/varlociraptor.yaml"
    params:
        event="PRESENT",
    shell:
        """
        PIPELINE_PATH=$(pwd)
        cd {input.varlo_path}
        cargo run --release -- filter-calls control-fdr --mode local-smart {input.calls} --events {params.event} --fdr 1 > {output}
        """


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
