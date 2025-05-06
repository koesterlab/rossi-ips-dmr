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
        git clone git@github.com:varlociraptor/varlociraptor.git varlociraptor
        cd varlociraptor
        git checkout methylation-paired-end
        """


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
        """ 
        samtools faidx {input}
        """


rule index_alignment:
    input:
        "resources/{platform}/{sample}.bam",
    output:
        "resources/{platform}/{sample}.bam.bai",
    conda:
        "../envs/samtools.yaml"
    threads: 10
    shell:
        """
        samtools index -@ {threads} {input}
        """


# Problem: The candidates span more than one chromosome... We would have to look at each chromosome indiviually
rule scatter_aligned_reads:
    input:
        alignment="resources/{platform}/{sample}.bam",
        candidate=lambda wildcards: "resources/candidates/candidates_{scatteritem}.bed",
    output:
        "resources/{platform}/{sample}/alignment_{scatteritem}.bam",
    log:
        "logs/scatter_aligned_reads_{platform}_{sample}_{scatteritem}.log",
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools view -b -L {input.candidate} {input.alignment} > {output}"
