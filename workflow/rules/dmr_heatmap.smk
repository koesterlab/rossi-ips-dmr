rule compute_chipseeker_regions:
    input:
        metilene="results/dmr_calls/{group2}/metilene_output.bed",
    output:
        chipseeker="results/dmr_calls/{group2}/chipseeker_regions.tsv",
    conda:
        "../envs/chipseeker.yaml"
    script:
        "../scripts/chipseeker.R"


rule create_heatmap:
    input:
        expand(
            "results/dmr_calls/{sample}-ref-sorted/chipseeker_regions.tsv",
            sample=["BC02", "BC03", "BC04"],
        ),
    output:
        "results/dmr_calls/heatmap.png",
    conda:
        "../envs/plot.yaml"
    script:
        "../scripts/heatmap.py"
