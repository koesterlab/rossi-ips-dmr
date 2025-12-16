# Needs a fasta with >chr1 instead of >1
rule modkit_compute_methylation:
    input:
        alignment="resources/nanopore/{sample}.bam",
        alignment_index="resources/nanopore/{sample}.bam.bai",
        genome="resources/genome.fasta",
    output:
        "results/nanopore/{caller}/meth_calling/{sample}/alignments_CpG.combined.bed",
    conda:
        "../envs/modkit.yaml"
    log:
        "logs/modkit_compute_methylation/{caller}_{sample}.log",
    resources:
        mem_mb=16000
    shell:
        """
        export PATH=$PATH:~/.cargo/bin 2> {log}
        export PATH=$PATH:/homes/aprinz/.cargo/bin 2> {log}
        modkit pileup {input.alignment} {output} --cpg --ref {input.genome} --force-allow-implicit --combine-strands --log-filepath {log} 2> {log}
        """


rule modkit_rename_output:
    input:
        "results/nanopore/{caller}/meth_calling/{sample}/alignments_CpG.combined.bed",
    output:
        "results/nanopore/{caller}/meth_calling/{sample}/modkit.vcf",
    log:
        "logs/modkit_rename_output/{sample}_{caller}.log",
    shell:
        "mv {input} {output} 2> {log}"
