sample_tsv_path = config["sample_path"]
chromosome_conf = config["resources"]["ref"]


def read_sample_tsv(sample_tsv_path):
    samples = {}
    with open(sample_tsv_path, "r") as file:
        next(file)
        for line in file:
            name, path, sequencer = line.strip().split("\t")
            samples[name] = (path, sequencer)
    return samples


samples = read_sample_tsv(sample_tsv_path)


def get_bioc_species_name():
    first_letter = config["resources"]["ref"]["species"][0]
    subspecies = config["resources"]["ref"]["species"].split("_")[1]
    return first_letter + subspecies


def all_input(wildcards):
    wanted_input = []
    # wanted_input.extend(
    #     expand(
    #         "results/{platform}/dmr_calls/heatmaps/{type}.html",
    #         platform=[platform for platform in config["platforms"]],
    #         type=[
    #             "distal_intergenic",
    #             "promoter",
    #             "intron",
    #             "exon",
    #             "3_utr",
    #             "5_utr",
    #             "downstream",
    #         ],
    #     )
    # )

    wanted_input.extend(
        expand(
            "results/{platform}/dmr_calls/{group2}/plots/dmr_qval.0.05.pdf",
            platform=[platform for platform in config["platforms"]],
            group2=[
                sample for sample in samples.keys() if sample != config["ref_sample"]
            ],
        ),
    )
    # wanted_input.extend(
    #     expand(
    #         "results/{platform}/dmr_calls/{group2}/plots/pvals.html",
    #         group2=[
    #             sample for sample in samples.keys() if sample != config["ref_sample"]
    #         ],
    #     ),
    # )
    wanted_input.extend(
        expand(
            "results/{platform}/datavzrd-report/{group2}",
            platform=[platform for platform in config["platforms"]],
            group2=[
                sample for sample in samples.keys() if sample != config["ref_sample"]
            ],
        ),
    )

    # Plots paper
    wanted_input.extend(
        expand(
            "results/{platform}/plots_paper/{group2}/scatter_plot.html",
            platform=[platform for platform in config["platforms"]],
            group2=[
                sample for sample in samples.keys() if sample != config["ref_sample"]
            ],
        ),
    )
    wanted_input.extend(
        expand(
            "results/{platform}/plots_paper/endo_meso/scatter_plot.html",
            platform=[platform for platform in config["platforms"]],
        )
    )
    wanted_input.extend(
        expand(
            "results/{platform}/plots_paper/pluripotency_score_all.html",
            platform=[platform for platform in config["platforms"]],
        )
    )
    wanted_input.extend(
        expand(
            "results/{platform}/plots_paper/pluripotency_score_psc.html",
            platform=[platform for platform in config["platforms"]],
        )
    )
    wanted_input.extend(
        expand(
            "results/comp_pb_np/meth_comp_pb_np_{group}.html",
            group=[sample for sample in samples.keys()],
        )
    )

    return wanted_input
