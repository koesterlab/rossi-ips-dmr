
rule get_fastq_se_gz:
    output:
        "resources/rna_seq/fastqs/{accession}.fastq.gz"
    log:
        "logs/get_fastq_se_gz/{accession}.gz.log"
    params:
        extra="--skip-technical"
    threads: 6
    wrapper:
        "v7.6.0/bio/sra-tools/fasterq-dump"

rule download_all_rnas:
    input:
        expand("resources/rna_seq/fastqs/{accession}.fastq.gz", accession=config["rna_accessions"])
    output:
        "rna_download_complete.txt"
    shell:
        """
        touch {output}
        """