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
    shell:
        "samtools faidx {input} 2> {log}"


# Indexes both the full alignments and their per-chunk subsets.
rule index_alignment:
    input:
        "resources/{platform}/{bam}.bam",
    output:
        "resources/{platform}/{bam}.bam.bai",
    conda:
        "../envs/samtools.yaml"
    threads: 10
    log:
        "logs/index_alignment/{platform}/{bam}.log",
    shell:
        "samtools index -@ {threads} {input} 2> {log}"



if not config["use_precomputed_calls"]:

    # Only keep reads overlapping the candidates of one chunk, so that each
    # varlociraptor job only has to process a small alignment file.
    rule scatter_aligned_reads:
        input:
            alignment="resources/{platform}/{germ_layer}.bam",
            candidate="resources/candidates/candidates_{scatteritem}.bed",
        output:
            "resources/{platform}/{germ_layer}/alignment_{scatteritem}.bam",
        log:
            "logs/scatter_aligned_reads/{platform}_{germ_layer}_{scatteritem}.log",
        conda:
            "../envs/samtools.yaml"
        shell:
            "samtools view -b -L {input.candidate} {input.alignment} > {output} 2> {log}"
