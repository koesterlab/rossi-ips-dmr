
scattergather:
    split_candidates=config["scatter_items"],


rule find_candidates:
    input:
        "resources/genome.fasta",
    output:
        "resources/candidates.bcf",
    log:
        "logs/find_candidates.log",
    conda:
        "../envs/varlociraptor.yaml"
    params:
        varlo_path=config["varlo_path"],
        pipeline_path=config["pipeline_path"],
    shell:
        """ 
        cd {params.varlo_path}
        cargo run -- methylation-candidates {params.pipeline_path}{input} {params.pipeline_path}{output}
        """


rule split_candidates:
    input:
        "resources/candidates.bcf",
    output:
        scatter.split_candidates("resources/candidates_{scatteritem}.bcf"),
    log:
        "logs/split_candidates.log",
    conda:
        "../envs/rbt.yaml"
    shell:
        "rbt vcf-split {input} {output}"


rule compute_meth_observations:
    input:
        genome="resources/genome.fasta",
        genomeIndex="resources/genome.fasta.fai",
        # alignments="resources/alignments/BC04-ref-sorted_debug.bam",
        alignments="resources/alignments/{sample}.bam",
        # alignments_index="resources/alignments/BC04-ref-sorted_debug.bam.bai",
        alignments_index="resources/alignments/{sample}.bam.bai",
        candidates="resources/candidates_{scatteritem}.bcf",
        # candidates="resources/candidates_debug_{scatteritem}.bcf",
    output:
        "results/{sample}/normal_{scatteritem}.bcf",
        # "results/normal_{scatteritem}.bcf",
    log:
        "logs/compute_meth_observations_{sample}_{scatteritem}.log",
    conda:
        "../envs/varlociraptor.yaml"
    params:
        varlo_path=config["varlo_path"],
        pipeline_path=config["pipeline_path"],
        sequencer=lambda wildcards: samples[wildcards.sample][1],
        # sequencer="Nanopore",
    shell:
        """ 
        cd {params.varlo_path}
        cargo run --release -- preprocess variants {params.pipeline_path}{input.genome} --candidates {params.pipeline_path}{input.candidates} --bam {params.pipeline_path}{input.alignments} --read-type {params.sequencer} > {params.pipeline_path}{output}
        """


rule call_methylation:
    input:
        preprocess_obs="results/{sample}/normal_{scatteritem}.bcf",
        scenario="resources/scenario.yaml",
    output:
        "results/{sample}/calls_{scatteritem}.bcf",
    log:
        "logs/call_methylation_{sample}_{scatteritem}.log",
    conda:
        "../envs/varlociraptor.yaml"
    params:
        varlo_path=config["varlo_path"],
        pipeline_path=config["pipeline_path"],
    shell:
        """ 
        cd {params.varlo_path}
        cargo run --release -- call variants --omit-strand-bias generic --scenario {params.pipeline_path}{input.scenario} --obs normal={params.pipeline_path}{input.preprocess_obs} > {params.pipeline_path}{output}
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
