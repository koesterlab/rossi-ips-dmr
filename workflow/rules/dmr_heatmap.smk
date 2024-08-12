config_ref = config["resources"]["ref"]


rule create_heatmap:
    input:
        expand(
            "results/dmr_calls/{sample}-ref-sorted/genes_transcripts/chipseeker_regions.tsv",
            sample=["BC02", "BC03", "BC04"],
        ),
    output:
        report(
            "results/dmr_calls/heatmaps/{type}.png",
            category="DMR plots",
            subcategory="heatmaps",
            labels=lambda wildcards: {
                "genetic element": wildcards.type,
            },
        ),
    conda:
        "../envs/plot.yaml"
    script:
        "../scripts/heatmap.py"
