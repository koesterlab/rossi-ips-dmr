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


rule aligned_downsampled_index:
    input:
        "resources/alignments/{sample}.bam",
    output:
        "resources/alignments/{sample}.bam.bai",
    log:
        "logs/aligned_downsampled_index_{sample}.log",
    conda:
        "../envs/samtools.yaml"
    threads: 10
    shell:
        """
        samtools index -@ {threads} {input}
        """
