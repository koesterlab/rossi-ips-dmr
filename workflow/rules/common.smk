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
            f"results/{platform}/{caller}/dmr_calls/heatmaps/{annotation}.png"
            for platform in config["meth_caller"].keys()
            for caller in config["meth_caller"].get(platform, [])
            for annotation in [
                "distal_intergenic",
                "promoter",
                "intron",
                "exon",
                "3_utr",
                "5_utr",
                # "downstream",
            ]
        ]
    )

    wanted_input.extend(
        [
            f"results/{platform}/{caller}/{rna_data}/diffexp_vs_dmrs_{tf}_tfs_{annotation_type}"
            for platform in config["meth_caller"].keys()
            for caller in config["meth_caller"].get(platform, [])
            for rna_data in ["rna_new", "rna_old"]
            for tf in ["no", "with"]
            for annotation_type in [
                # "distal_intergenic",
                "promoter",
                # "intron",
                # "exon",
                # "3_utr",
                # "5_utr",
                # "downstream",
            ]
        ]
    )

    wanted_input.extend(
        [
            f"results/{platform}/{caller}/{rna_data}/pathways/{germ_layer}-gene_set_{annotation_type}-{func}"
            for platform in config["meth_caller"].keys()
            for caller in config["meth_caller"].get(platform, [])
            for rna_data in ["rna_new", "rna_old"]
            for germ_layer in ["endoderm", "mesoderm", "ectoderm", "all"]
            for annotation_type in [
                # "distal_intergenic",
                "promoter",
                # "intron",
                # "exon",
                # "3_utr",
                # "5_utr",
                # "downstream",
            ]
            for func in ["mf", "bp", "cc", "go"]
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
            f"results/{platform}/{caller}/plots_paper/{group2}/scatter_plot.png"
            for platform in config["meth_caller"].keys()
            for caller in config["meth_caller"].get(platform, [])
            for group2 in [s for s in samples.keys() if s != config["ref_sample"]]
        ]
    )

    wanted_input.extend(
        [
            f"results/{platform}/{caller}/plots_paper/endo_meso/scatter_plot.png"
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
    return wanted_input
