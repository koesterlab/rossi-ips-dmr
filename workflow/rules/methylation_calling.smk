
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
        PIPELINE_PATH=$(pwd)
        cd {input.varlo_path}
        cargo run -- methylation-candidates $PIPELINE_PATH/{input.fasta} $PIPELINE_PATH/{output}
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


rule compute_meth_observations:
    input:
        varlo_path="resources/tools/varlociraptor",
        genome="resources/genome.fasta",
        genomeIndex="resources/genome.fasta.fai",
        alignments="resources/alignments/{sample}.bam",
        alignments_index="resources/alignments/{sample}.bam.bai",
        candidates="resources/candidates/candidates_{scatteritem}.bcf",
    output:
        "results/{sample}/normal_{scatteritem}.bcf",
    log:
        "logs/compute_meth_observations_{sample}_{scatteritem}.log",
    conda:
        "../envs/varlociraptor.yaml"
    params:
        sequencer=lambda wildcards: samples[wildcards.sample][1],
    shell:
        """ 
        PIPELINE_PATH=$(pwd)
        cd {input.varlo_path}
        cargo run --release -- preprocess variants $PIPELINE_PATH/{input.genome} --candidates $PIPELINE_PATH/{input.candidates} --bam $PIPELINE_PATH/{input.alignments} --read-type {params.sequencer} --max-depth 1000 > $PIPELINE_PATH/{output}
        """


rule call_methylation:
    input:
        varlo_path="resources/tools/varlociraptor",
        preprocess_obs="results/{sample}/normal_{scatteritem}.bcf",
        scenario="resources/scenario.yaml",
    output:
        "results/{sample}/calls_{scatteritem}.bcf",
    log:
        "logs/call_methylation_{sample}_{scatteritem}.log",
    conda:
        "../envs/varlociraptor.yaml"
    shell:
        """ 
        PIPELINE_PATH=$(pwd)
        cd {input.varlo_path}
        cargo run --release -- call variants --omit-strand-bias generic --scenario $PIPELINE_PATH/{input.scenario} --obs normal=$PIPELINE_PATH/{input.preprocess_obs} > $PIPELINE_PATH/{output}
        """


rule calls_to_vcf:
    input:
        "results/{sample}/calls_{scatteritem}.bcf",
    output:
        "results/{sample}/calls_{scatteritem}.vcf",
    conda:
        "../envs/samtools.yaml"
    log:
        "logs/convert_to_vcf_{sample}_{scatteritem}.log",
    threads: 10
    shell:
        """
        bcftools view --threads {threads} {input} -o {output}
        """


rule gather_calls:
    input:
        gather.split_candidates("results/{{sample}}/calls_{scatteritem}.vcf"),
    output:
        "results/{sample}/calls.vcf",
    log:
        "logs/gather_calls_{sample}.log",
    conda:
        "../envs/cat.yaml"
    shell:
        "cat {input} > {output}"
