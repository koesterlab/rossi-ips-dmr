chromosome_conf = config["sample"]


rule get_genome:
    output:
        "resources/genome.fasta",
    params:
        species=chromosome_conf["species"],
        datatype=chromosome_conf["datatype"],
        build=chromosome_conf["build"],
        release=chromosome_conf["release"],
    log:
        "logs/get_genome.log",
    cache: "omit-software"  # save space and time with between workflow caching (see docs)
    wrapper:
        "v2.3.2/bio/reference/ensembl-sequence"


rule genome_index:
    input:
        "resources/genome.fasta",
    output:
        "resources/genome.fasta.fai",
    log:
        "logs/genome_index.log",
    conda:
        "../envs/samtools.yaml"
    params:
        pipeline_path=config["pipeline_path"],
    shell:
        """ 
        samtools faidx {params.pipeline_path}{input}
        """


rule aligned_downsampled_index:
    input:
        "resources/alignments/{sample}.bam",
    output:
        "resources/alignments/{sample}.bam.bai",
    log:
        "logs/aligned_downsampled_index_{sample}.log",
    conda:
        "../envs/samtools.yaml"
    params:
        pipeline_path=config["pipeline_path"],
    threads: 10
    shell:
        """
        samtools index -@ {threads} {params.pipeline_path}{input}
        """
