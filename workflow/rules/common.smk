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

    wanted_input.extend(
        [
            f"results/{platform}/{caller}/dmr_calls/heatmaps/{type}.png"
            for platform in config["meth_caller"].keys()
            for caller in config["meth_caller"].get(platform, [])
            for type in [
                "distal_intergenic",
                "promoter",
                "intron",
                "exon",
                "3_utr",
                "5_utr",
                "downstream",
            ]
        ]
    )

    wanted_input.extend(
        [
            f"results/{platform}/{caller}/dmr_calls/{group2}/plots/dmr_qval.0.05.pdf"
            for platform in config["meth_caller"].keys()
            for caller in config["meth_caller"].get(platform, [])
            for group2 in [s for s in samples.keys() if s != config["ref_sample"]]
        ]
    )

    wanted_input.extend(
        [
            f"results/{platform}/{caller}/datavzrd-report/{group2}"
            for platform in config["meth_caller"].keys()
            for caller in config["meth_caller"].get(platform, [])
            for group2 in [s for s in samples.keys() if s != config["ref_sample"]]
        ]
    )

    wanted_input.extend(
        [
            f"results/{platform}/{caller}/rna_seq/datavzrd/{group2}"
            for platform in config["meth_caller"].keys()
            for caller in config["meth_caller"].get(platform, [])
            for group2 in [s for s in samples.keys() if s != config["ref_sample"]]
        ]
    )

    wanted_input.extend(
        [
            f"results/{platform}/{caller}/plots_paper/{group2}/scatter_plot.html"
            for platform in config["meth_caller"].keys()
            for caller in config["meth_caller"].get(platform, [])
            for group2 in [s for s in samples.keys() if s != config["ref_sample"]]
        ]
    )

    wanted_input.extend(
        [
            f"results/{platform}/{caller}/plots_paper/endo_meso/scatter_plot.html"
            for platform in config["meth_caller"].keys()
            for caller in config["meth_caller"].get(platform, [])
        ]
    )

    wanted_input.extend(
        [
            f"results/{platform}/{caller}/plots_paper/pluripotency_score_all.html"
            for platform in config["meth_caller"].keys()
            for caller in config["meth_caller"].get(platform, [])
        ]
    )

    wanted_input.extend(
        [
            f"results/{platform}/{caller}/plots_paper/pluripotency_score_psc.html"
            for platform in config["meth_caller"].keys()
            for caller in config["meth_caller"].get(platform, [])
            for group2 in [s for s in samples.keys() if s != config["ref_sample"]]
        ]
    )

    # wanted_input.extend([
    #         f"results/comp_pb_np/meth_comp_pb_np_{group}.png"
    #         for group in [s for s in samples.keys() ]
    #     ])

    return wanted_input
