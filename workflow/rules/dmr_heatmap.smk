rule create_heatmap:
    input:
        expand(
            "results/{{platform}}/dmr_calls/{sample}/genes_transcripts/chipseeker.tsv",
            sample=[
                sample for sample in samples.keys() if sample != config["ref_sample"]
            ],
        ),
    output:
        report(
            "results/{platform}/dmr_calls/heatmaps/{type}.html",
            caption="../report/heatmap.rst",
            category="DMR plots",
            subcategory="Heatmaps",
            labels=lambda wildcards: {
                "platform": wildcards.platform,
                "genetic element": wildcards.type,
            },
        ),
    conda:
        "../envs/plot.yaml"
    script:
        "../scripts/heatmap.py"
