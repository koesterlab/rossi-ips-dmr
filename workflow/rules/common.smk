sample_tsv_path = config["sample_path"]
chromosome_conf = config["resources"]["ref"]

ALL_GERM_LAYERS = ["psc", "endoderm", "mesoderm", "ectoderm"]


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


def get_non_base_layers(base):
    """Return the 3 germ layers that are NOT the base experiment."""
    return [layer for layer in ALL_GERM_LAYERS if layer != base]


def get_base_experiments():
    """All 4 germ layers serve as base experiment once."""
    return ALL_GERM_LAYERS


def get_rna_data_values():
    """
    Return the two rna_data wildcard values.
    The {rna_data} wildcard is always strictly 'rna_old' or 'rna_new'.
    The base_level is carried by the separate {base} wildcard.
    The kallisto-sleuth module prefix is rna_{rna_data}/base_{base},
    so the diffexp tables live at:
      rna_old/base_{base}/results/tables/diffexp/...
      rna_new/base_{base}/results/tables/diffexp/...
    """
    return config["rna_data"]


def all_input(wildcards):
    wanted_input = []

    # DMR heatmaps – one set per base experiment
    wanted_input.extend(
        [
            f"results/{platform}/{caller}/base_{base}/dmr_calls/heatmaps/{annotation}.png"
            for platform in config["meth_caller"].keys()
            for caller in config["meth_caller"].get(platform, [])
            for base in get_base_experiments()
            for annotation in [
                "distal_intergenic",
                "promoter",
                "intron",
                "exon",
                "3_utr",
                "5_utr",
            ]
        ]
    )

    # DMR vs DiffExp comparisons – 4 bases × 2 RNA datasets × 2 tf modes × annotation types
    wanted_input.extend(
        [
            f"results/{platform}/{caller}/base_{base}/{rna_data}/diffexp_vs_dmrs_{tf}_tfs_{annotation_type}"
            for platform in config["meth_caller"].keys()
            for caller in config["meth_caller"].get(platform, [])
            for base in get_base_experiments()
            for rna_data in get_rna_data_values()
            # for tf in ["no"]
            for tf in ["no", "with"]
            for annotation_type in [
                "promoter",
            ]
        ]
    )

    # fgsea pathway enrichment – 4 bases × 2 RNA datasets × non-base layers + "all" × annotation × func
    wanted_input.extend(
        [
            f"results/{platform}/{caller}/base_{base}/{rna_data}/pathways/{germ_layer}-gene_set_{annotation_type}-{func}"
            for platform in config["meth_caller"].keys()
            for caller in config["meth_caller"].get(platform, [])
            for base in get_base_experiments()
            for rna_data in get_rna_data_values()
            for germ_layer in get_non_base_layers(base) + ["all"]
            for annotation_type in [
                "promoter",
            ]
            for func in ["mf", "bp", "cc", "go"]
        ]
    )

    # Metilene plots (PDF) – one per base × non-base group
    wanted_input.extend(
        [
            f"results/{platform}/{caller}/base_{base}/dmr_calls/{group2}/plots/dmr_qval.0.05.pdf"
            for platform in config["meth_caller"].keys()
            for caller in config["meth_caller"].get(platform, [])
            for base in get_base_experiments()
            for group2 in get_non_base_layers(base)
        ]
    )

    # datavzrd annotation reports – one per base × non-base group
    wanted_input.extend(
        [
            f"results/{platform}/{caller}/base_{base}/datavzrd-report/{group2}"
            for platform in config["meth_caller"].keys()
            for caller in config["meth_caller"].get(platform, [])
            for base in get_base_experiments()
            for group2 in get_non_base_layers(base)
        ]
    )

    # Scatter plots – one per base × non-base group
    wanted_input.extend(
        [
            f"results/{platform}/{caller}/base_{base}/plots_paper/{group2}/scatter_plot.png"
            for platform in config["meth_caller"].keys()
            for caller in config["meth_caller"].get(platform, [])
            for base in get_base_experiments()
            for group2 in get_non_base_layers(base)
        ]
    )

    # Endo-meso scatter plots – only relevant when neither endoderm nor mesoderm is the base
    wanted_input.extend(
        [
            f"results/{platform}/{caller}/plots_paper/endo_meso/scatter_plot.png"
            for platform in config["meth_caller"].keys()
            for caller in config["meth_caller"].get(platform, [])
            for base in get_base_experiments()
            if base not in ("endoderm", "mesoderm")
        ]
    )

    # Pluripotency score heatmaps – one per base (data source is the same parquet,
    # kept scoped per base for consistency)
    wanted_input.extend(
        [
            f"results/{platform}/{caller}/plots_paper/pluripotency_score_all.html"
            for platform in config["meth_caller"].keys()
            for caller in config["meth_caller"].get(platform, [])
            for base in get_base_experiments()
        ]
    )

    return wanted_input
