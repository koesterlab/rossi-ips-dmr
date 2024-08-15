rule create_heatmap:
    input:
        expand(
            "results/dmr_calls/{sample}/genes_transcripts/chipseeker.tsv",
            sample=[
                sample for sample in samples.keys() if sample != config["ref_sample"]
            ],
        ),
    output:
        report(
            "results/dmr_calls/heatmaps/{type}.png",
            caption="../report/heatmap.rst",
            category="DMR plots",
            subcategory="Heatmaps",
            labels=lambda wildcards: {
                "genetic element": wildcards.type,
            },
        ),
    conda:
        "../envs/plot.yaml"
    script:
        "../scripts/heatmap.py"
