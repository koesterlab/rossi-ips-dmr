import yaml


def get_bioc_species_pkg():
    """Get the bioconductor annotation package name for the species in config.yaml"""
    species_letters = get_bioc_species_name()[0:2].capitalize()
    return "org.{species}.eg.db".format(species=species_letters)


def render_enrichment_env():
    species_pkg = f"bioconductor-{get_bioc_species_pkg()}"
    with open(workflow.source_path("../envs/enrichment.yaml")) as f:
        env = yaml.load(f, Loader=yaml.SafeLoader)
    env["dependencies"].append(species_pkg)
    env_path = Path("resources/envs/enrichment.yaml")
    rendered = yaml.dump(env)
    # Only write if changed, so that the file is not touched on every parse
    if not env_path.exists() or env_path.read_text() != rendered:
        env_path.parent.mkdir(parents=True, exist_ok=True)
        env_path.write_text(rendered)
    return env_path.absolute()


func_to_names = {
    "mf": "molecular_function",
    "bp": "biological_process",
    "cc": "cellular_component",
    "go": "all",
}

bioc_species_pkg = get_bioc_species_pkg()
enrichment_env = render_enrichment_env()

# fgsea is run per non-base germ layer and for all layers together
FGSEA_GERM_LAYERS = "|".join(ALL_GERM_LAYERS + ["all"])


rule fgsea_dmr_vs_diffexp:
    input:
        # The TF-adjusted table provides the ranking column (ranked_meth_diffexp)
        diffexp_vs_dmrs_promoter="results/{platform}/{caller}/base_{base}/{rna_data}/tfs/diffexp_vs_dmrs_{annotation_type}.tsv",
        gene_sets=lambda wildcards: config["fgsea"][f"gene_sets_{wildcards.func}"],
        common_src=workflow.source_path("../scripts/common.R"),
    output:
        enrichment="results/{platform}/{caller}/base_{base}/{rna_data}/rna_seq_comp/{germ_layer}-all-gene-sets-{annotation_type}-{func}.tsv",
        rank_ties="results/{platform}/{caller}/base_{base}/{rna_data}/rna_seq_comp/{germ_layer}-rank-ties-{annotation_type}-{func}.tsv",
        significant="results/{platform}/{caller}/base_{base}/{rna_data}/rna_seq_comp/{germ_layer}-sig-gene-sets-{annotation_type}-{func}.tsv",
        plot="results/{platform}/{caller}/base_{base}/{rna_data}/rna_seq_comp/{germ_layer}-table-plot-{annotation_type}-{func}.pdf",
        plot_collapsed="results/{platform}/{caller}/base_{base}/{rna_data}/rna_seq_comp/{germ_layer}-collapsed_pathways.table-plot-{annotation_type}-{func}.pdf",
    wildcard_constraints:
        germ_layer=FGSEA_GERM_LAYERS,
    params:
        bioc_species_pkg=bioc_species_pkg,
        gene_set_fdr=config["fgsea"]["fdr_gene_set"],
        eps=config["fgsea"]["eps"],
        germ_layer=lambda wildcards: wildcards.germ_layer,
    conda:
        enrichment_env
    log:
        "logs/fgsea_dmr_vs_diffexp/{platform}_{caller}_{base}_{rna_data}_{germ_layer}_{annotation_type}_{func}.log",
    threads: 25
    script:
        "../scripts/fgsea.R"


rule fgsea_datavzrd:
    input:
        config=workflow.source_path("../resources/fgsea.yaml"),
        comp="results/{platform}/{caller}/base_{base}/{rna_data}/rna_seq_comp/{germ_layer}-all-gene-sets-{annotation_type}-{func}.tsv",
    output:
        report(
            directory(
                "results/{platform}/{caller}/base_{base}/{rna_data}/pathways/{germ_layer}-gene_set_{annotation_type}-{func}"
            ),
            caption="../report/diffexp_vs_dmrs.rst",
            htmlindex="index.html",
            category=lambda wildcards: f"DiffExp-DMRs Pathways - {wildcards.rna_data}",
            subcategory=lambda wildcards: f"Base: {wildcards.base}",
            labels=lambda wildcards: {
                "comparison": wildcards.germ_layer,
                "func": func_to_names[wildcards.func],
            },
        ),
    wildcard_constraints:
        germ_layer=FGSEA_GERM_LAYERS,
    log:
        "logs/fgsea_datavzrd/{platform}_{caller}_{base}_{rna_data}_{germ_layer}_{annotation_type}_{func}.log",
    wrapper:
        "v9.4.1/utils/datavzrd"
