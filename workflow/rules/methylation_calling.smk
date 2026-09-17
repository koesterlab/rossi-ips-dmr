scattergather:
    split_candidates=config["scatter_items"],


if not config["use_precomputed_calls"]:

    rule find_candidates:
        input:
            fasta="resources/genome.fasta",
        output:
            "resources/candidates/candidates.bcf",
        log:
            "logs/varlociraptor/find_candidates.log",
        conda:
            "../envs/varlociraptor.yaml"
        shell:
            "varlociraptor methylation-candidates {input.fasta} {output} --motifs CG 2> {log}"

    rule split_candidates:
        input:
            "resources/candidates/candidates.bcf",
        output:
            scatter.split_candidates("resources/candidates/candidates_{scatteritem}.bcf"),
        log:
            "logs/varlociraptor/split_candidates.log",
        conda:
            "../envs/rbt.yaml"
        shell:
            "rbt vcf-split {input} {output} 2> {log}"

    rule compute_meth_observations:
        input:
            genome="resources/genome.fasta",
            genome_index="resources/genome.fasta.fai",
            alignments="resources/{platform}/{germ_layer}/alignment_{scatteritem}.bam",
            alignments_index="resources/{platform}/{germ_layer}/alignment_{scatteritem}.bam.bai",
            candidates="resources/candidates/candidates_{scatteritem}.bcf",
        output:
            "results/{platform}/varlo/meth_calling/{germ_layer}/normal_{scatteritem}.bcf",
        conda:
            "../envs/varlociraptor.yaml"
        log:
            "logs/varlociraptor/compute_meth_observations/{platform}_{germ_layer}_{scatteritem}.log",
        shell:
            "varlociraptor preprocess variants {input.genome} --candidates {input.candidates} "
            "--bam {input.alignments} --methylation-read-type annotated --max-depth 1000 "
            "> {output} 2> {log}"

    rule call_methylation_single:
        input:
            preprocess_obs="results/{platform}/varlo/meth_calling/{germ_layer}/normal_{scatteritem}.bcf",
            scenario=workflow.source_path("../resources/scenarios/scenario.yaml"),
        output:
            "results/{platform}/varlo/meth_calling/{germ_layer}/calls_{scatteritem}.bcf",
        conda:
            "../envs/varlociraptor.yaml"
        wildcard_constraints:
            platform="pacbio|nanopore",
        log:
            "logs/call_methylation_single/{platform}_{germ_layer}_{scatteritem}.log",
        shell:
            "varlociraptor call variants generic --scenario {input.scenario} "
            "--obs normal={input.preprocess_obs} > {output} 2> {log}"

    # Joint calling of PacBio and Nanopore observations.
    rule call_methylation:
        input:
            pb="results/pacbio/varlo/meth_calling/{germ_layer}/normal_{scatteritem}.bcf",
            np="results/nanopore/varlo/meth_calling/{germ_layer}/normal_{scatteritem}.bcf",
            scenario=workflow.source_path("../resources/scenarios/scenario_common.yaml"),
        output:
            "results/platforms_combined/varlo/meth_calling/{germ_layer}/calls_{scatteritem}.bcf",
        conda:
            "../envs/varlociraptor.yaml"
        log:
            "logs/varlociraptor/compute_meth_together/{germ_layer}_{scatteritem}.log",
        resources:
            mem_mb=128000,
        shell:
            "varlociraptor call variants generic --scenario {input.scenario} "
            "--obs pacbio={input.pb} nanopore={input.np} > {output} 2> {log}"

    # Control the FDR separately for PRESENT and ABSENT calls; both sets are
    # merged again afterwards.
    rule filter_calls:
        input:
            "results/{platform}/varlo/meth_calling/{germ_layer}/calls_{scatteritem}.bcf",
        output:
            "results/{platform}/varlo/meth_calling/{germ_layer}/{fdr}/calls_{scatteritem}.filtered_{event}.bcf",
        wildcard_constraints:
            event="absent|present",
        conda:
            "../envs/varlociraptor.yaml"
        params:
            event=lambda wildcards: wildcards.event.upper(),
        log:
            "logs/varlociraptor/filter_calls/{platform}_{germ_layer}_{scatteritem}_{fdr}_{event}.log",
        shell:
            "varlociraptor filter-calls control-fdr --mode local-smart {input} "
            "--events {params.event} --fdr {wildcards.fdr} > {output} 2> {log}"

    rule concatenate_filtered_calls:
        input:
            absent="results/{platform}/varlo/meth_calling/{germ_layer}/{fdr}/calls_{scatteritem}.filtered_absent.bcf",
            absent_index="results/{platform}/varlo/meth_calling/{germ_layer}/{fdr}/calls_{scatteritem}.filtered_absent.bcf.csi",
            present="results/{platform}/varlo/meth_calling/{germ_layer}/{fdr}/calls_{scatteritem}.filtered_present.bcf",
            present_index="results/{platform}/varlo/meth_calling/{germ_layer}/{fdr}/calls_{scatteritem}.filtered_present.bcf.csi",
        output:
            "results/{platform}/varlo/meth_calling/{germ_layer}/{fdr}/calls_{scatteritem}.filtered.bcf",
        conda:
            "../envs/samtools.yaml"
        log:
            "logs/varlociraptor/concatenate_filtered_calls/{platform}_{germ_layer}_{scatteritem}_{fdr}.log",
        shell:
            # Remove records that occur in both files
            "(bcftools concat -a {input.absent} {input.present} | "
            "bcftools norm -d exact -o {output}) 2> {log}"

    rule gather_calls:
        input:
            gather.split_candidates(
                "results/{{platform}}/varlo/meth_calling/{{germ_layer}}/{{fdr}}/calls_{scatteritem}.filtered.bcf"
            ),
        output:
            "results/{platform}/varlo/meth_calling/{germ_layer}/varlo_{fdr}.bcf",
        conda:
            "../envs/samtools.yaml"
        log:
            "logs/varlociraptor/gather_calls/{platform}_{germ_layer}_{fdr}.log",
        shell:
            "bcftools concat {input} -o {output} 2> {log}"


rule index_bcf:
    input:
        "{bcf}.bcf",
    output:
        "{bcf}.bcf.csi",
    log:
        "logs/index_bcf/{bcf}.log",
    conda:
        "../envs/samtools.yaml"
    shell:
        "bcftools index -c {input} 2> {log}"


rule prepare_wsabi:
    input:
        "results/{platform}/varlo/meth_calling/{germ_layer}/varlo_0.05.bcf",
    output:
        "results/wsabi/{platform}_{germ_layer}.tsv.gz",
    conda:
        "../envs/pysam.yaml"
    log:
        "logs/varlociraptor/prepare_wsabi/{platform}_{germ_layer}.log",
    script:
        "../scripts/bcf_to_wsabi_tsv.py"


rule df_from_calls:
    input:
        # One call file (plus index) per germ layer, named by layer, e.g. "psc", "psc_index"
        **{
            f"{layer}{key_suffix}": f"results/{{platform}}/{{caller}}/meth_calling/{layer}/{{caller}}_{{fdr}}.bcf{file_suffix}"
            for layer in ALL_GERM_LAYERS
            for key_suffix, file_suffix in [("", ""), ("_index", ".csi")]
        },
    output:
        "results/{platform}/{caller}/meth_calling/calls_{fdr}.parquet",
    conda:
        "../envs/pysam.yaml"
    log:
        "logs/varlociraptor/df_from_calls/{platform}_{caller}_{fdr}.log",
    resources:
        mem_mb=32000,
    script:
        "../scripts/df_from_calls.py"


rule annotate_methylation:
    input:
        methylation="results/{platform}/{caller}/meth_calling/calls_{fdr}.parquet",
        txdb="resources/ref/txdb.db",
        txnames="resources/ref/txnames.rds",
    output:
        chipseeker="results/{platform}/{caller}/meth_calling/chipseeker_{fdr}.tsv",
    log:
        "logs/annotate_methylation/{platform}_{caller}_{fdr}.log",
    conda:
        "../envs/chipseeker.yaml"
    script:
        "../scripts/chipseeker_methylation.R"
