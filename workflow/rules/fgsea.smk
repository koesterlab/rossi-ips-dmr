import yaml


def get_bioc_species_name():
    first_letter = config["resources"]["ref"]["species"][0]
    subspecies = config["resources"]["ref"]["species"].split("_")[1]
    return first_letter + subspecies


def get_bioc_species_pkg():
    """Get the package bioconductor package name for the the species in config.yaml"""
    species_letters = get_bioc_species_name()[0:2].capitalize()
    return "org.{species}.eg.db".format(species=species_letters)


def render_enrichment_env():
    species_pkg = f"bioconductor-{get_bioc_species_pkg()}"
    with open(workflow.source_path("../envs/enrichment.yaml")) as f:
        env = yaml.load(f, Loader=yaml.SafeLoader)
    env["dependencies"].append(species_pkg)
    env_path = Path("resources/envs/enrichment.yaml")
    env_path.parent.mkdir(parents=True, exist_ok=True)
    with open(env_path, "w") as f:
        yaml.dump(env, f)
    return env_path.absolute()


func_to_names = {
    "mf": "molecular_function",
    "bp": "biological_process",
    "cc": "cellular_component",
    "go": "all",
}

bioc_species_pkg = get_bioc_species_pkg()
enrichment_env = render_enrichment_env()


rule fgsea_dmr_vs_diffexp:
    input:
        # samples="results/sleuth/{model}.samples.tsv",
        # diffexp="results/tables/diffexp/{model}.genes-representative.diffexp.tsv",
        diffexp_vs_dmrs_promoter="results/{platform}/{caller}/{rna_data}/diffexp_vs_dmrs_{annotation_type}.tsv",
        gene_sets=lambda wildcards: config["fgsea"][f"gene_sets_{wildcards.func}"],
        common_src=workflow.source_path("../scripts/common.R"),
    output:
        enrichment="results/{platform}/{caller}/{rna_data}/{germ_layer}-all-gene-sets-{annotation_type}-{func}.tsv",
        rank_ties="results/{platform}/{caller}/{rna_data}/{germ_layer}-rank-ties-{annotation_type}-{func}.tsv",
        significant="results/{platform}/{caller}/{rna_data}/{germ_layer}-sig-gene-sets-{annotation_type}-{func}.tsv",
        plot="results/{platform}/{caller}/{rna_data}/{germ_layer}-table-plot-{annotation_type}-{func}.pdf",
        plot_collapsed="results/{platform}/{caller}/{rna_data}/{germ_layer}-collapsed_pathways.table-plot-{annotation_type}-{func}.pdf",
    params:
        bioc_species_pkg=bioc_species_pkg,
        # model=get_model,
        gene_set_fdr=config["fgsea"]["fdr_gene_set"],
        eps=config["fgsea"]["eps"],
        germ_layer=lambda wildcards: wildcards.germ_layer,
        # covariate=lambda w: config["diffexp"]["models"][w.model]["primary_variable"],
    conda:
        enrichment_env
    log:
        "logs/fgsea_dmr_vs_diffexp/{platform}_{caller}_{rna_data}_{germ_layer}_{annotation_type}_{func}.log",
    threads: 25
    script:
        "../scripts/fgsea.R"


rule fgsea_datavzrd:
    input:
        config=workflow.source_path("../resources/fgsea.yaml"),
        comp="results/{platform}/{caller}/{rna_data}/{germ_layer}-all-gene-sets-{annotation_type}-{func}.tsv",
    output:
        report(
            directory(
                "results/{platform}/{caller}/{rna_data}/pathways/{germ_layer}-gene_set_{annotation_type}-{func}"
            ),
            caption="../report/diffexp_vs_dmrs.rst",
            htmlindex="index.html",
            category="DiffExp-DMRs Comparison",
            subcategory=lambda wildcards: "pathways no transcription factors",
            labels=lambda wildcards: {
                "layer": wildcards.germ_layer,
                "func": wildcards.func,
            },
        ),
    log:
        "logs/diffexp_dmvzrd/diffexp_dmr_datavzrd/{platform}_{caller}_{rna_data}_{germ_layer}_{annotation_type}_{func}.log",
    wrapper:
        "641c90c4da86d4acf2022f347f3c8017334c0f44/utils/datavzrd"
