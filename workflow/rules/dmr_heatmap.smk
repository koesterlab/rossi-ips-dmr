rule create_heatmap:
    input:
        expand(
            "results/{{platform}}/{{caller}}/dmr_calls/{sample}/genes_transcripts/chipseeker_filtered.tsv",
            sample=[
                sample for sample in samples.keys() if sample != config["ref_sample"]
            ],
        ),
    output:
        report(
            "results/{platform}/{caller}/dmr_calls/heatmaps/{type}.png",
            caption="../report/heatmap.rst",
            category="DMR plots",
            subcategory=lambda wildcards: f"Heatmaps: {wildcards.platform} - {wildcards.caller}",
            labels=lambda wildcards: {
                "genetic element": wildcards.type,
            },
        ),
    conda:
        "../envs/plot.yaml"
    log:
        "logs/heatmap/{platform}/{caller}/{type}.log",
    script:
        "../scripts/heatmap.py"
