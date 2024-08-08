def all_input(wildcards):
    wanted_input = []
    wanted_input.extend(
        expand(
            "results/dmr_calls/heatmaps/{type}.png",
            type=[
                "distal_intergenic",
                "promoter",
                "intron",
                "exon",
                "3_utr",
                "5_utr",
                "downstream",
            ],
        )
    )
    wanted_input.extend(
        expand("results/{sample}/calls.vcf", sample=list(samples.keys())),
    )
    wanted_input.extend(
        expand(
            "results/dmr_calls/{group2}/plots/dmr_qval.0.05.pdf",
            group2=[
                sample for sample in samples.keys() if sample != config["ref_sample"]
            ],
        ),
    )
    wanted_input.extend(
        expand(
            "results/dmr_calls/{group2}/plots/pvals.png",
            group2=[
                sample for sample in samples.keys() if sample != config["ref_sample"]
            ],
        ),
    )
    wanted_input.extend(
        expand(
            "results/datavzrd-report/{group2}",
            group2=[
                sample for sample in samples.keys() if sample != config["ref_sample"]
            ],
        ),
    )
    return wanted_input
