chromosome_conf = config["resources"]["ref"]

ALL_GERM_LAYERS = ["psc", "endoderm", "mesoderm", "ectoderm"]
ACTIVE_BASE_COMP = ["psc"]
ANNOTATION_TYPES = ["distal_intergenic", "promoter", "intron", "exon", "3_utr", "5_utr", "unfiltered"]
# Use ["no", "with"] to additionally create the TF-adjusted comparison reports.
TF_MODES = ["no"]


# Constraining wildcards keeps the DAG construction fast: without them, generic
# patterns (e.g. "{bcf}.bcf.csi") match nearly every path and Snakemake has to
# try many candidate rules per file. Only names that are not used by the
# kallisto-sleuth module are constrained globally.
# wildcard_constraints:
#     platform="pacbio|nanopore|platforms_combined",
#     caller="|".join({c for callers in config["meth_caller"].values() for c in callers}),
#     base="|".join(ALL_GERM_LAYERS),
#     group2="|".join(ALL_GERM_LAYERS),
#     germ_layer="|".join(ALL_GERM_LAYERS),
#     fdr=r"\d+(?:\.\d+)?",
#     scatteritem=r"\d+-of-\d+",
#     rna_data="rna_old|rna_new",
#     annotation_type="|".join(ANNOTATION_TYPES),
#     plot_type="pdf|png|svg|html",


def get_bioc_species_name():
    first_letter = chromosome_conf["species"][0]
    subspecies = chromosome_conf["species"].split("_")[1]
    return first_letter + subspecies


def get_non_base_layers(base):
    """Return the 3 germ layers that are NOT the base experiment."""
    return [layer for layer in ALL_GERM_LAYERS if layer != base]


def platform_callers():
    """Yield all configured (platform, methylation caller) combinations."""
    for platform, callers in config["meth_caller"].items():
        for caller in callers:
            yield platform, caller


def chipseeker_tables(wildcards, fdr="0.05"):
    """Postprocessed ChIPseeker DMR annotations of all non-base layers."""
    return expand(
        "results/{platform}/{caller}/base_{base}/dmr_calls/{group2}/genes_transcripts/{fdr}/chipseeker_postprocessed.tsv",
        platform=wildcards.platform,
        caller=wildcards.caller,
        base=wildcards.base,
        group2=get_non_base_layers(wildcards.base),
        fdr=fdr,
    )


def all_input(wildcards):
    wanted_input = []

    for platform, caller in platform_callers():
        prefix = f"results/{platform}/{caller}"

        for base in ACTIVE_BASE_COMP:
            wanted_input += [
                f"{prefix}/base_{base}/{rna_data}/diffexp_vs_dmrs_{tf}_tfs_{annotation_type}"
                for rna_data in config["rna_data"]
                for tf in TF_MODES
                for annotation_type in ANNOTATION_TYPES
            ]
            for group2 in get_non_base_layers(base):
                wanted_input += [
                    f"{prefix}/base_{base}/dmr_calls/datavzrd-report/{group2}",
                ]
                if config["fgsea"]["activate"]:
                    wanted_input += [
                        f"{prefix}/base_{base}/{rna_data}/pathways/{germ_layer}-gene_set_promoter-{func}"
                        for rna_data in config["rna_data"]
                        for germ_layer in [group2, "all"]
                        for func in ["mf", "bp", "cc", "go"]
                    ]

        wanted_input += [
            f"{prefix}/plots_paper/pluripotency_score_all.html",
        ]

    wanted_input += [
        f"results/wsabi/{platform}_{layer}.tsv.gz"
        for platform in config["meth_caller"]
        for layer in ALL_GERM_LAYERS
    ]

    wanted_input += [
        "results/platforms_combined/varlo/plots_paper/scatter_comparison.pdf"
    ]

    # Duplicates (e.g. the "all" fgsea target) are harmless but slow down the DAG.
    return list(dict.fromkeys(wanted_input))
