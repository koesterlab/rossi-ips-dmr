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


bioc_species_pkg = get_bioc_species_pkg()
enrichment_env = render_enrichment_env()


rule fgsea_dmr_vs_diffexp:
    input:
        # samples="results/sleuth/{model}.samples.tsv",
        # diffexp="results/tables/diffexp/{model}.genes-representative.diffexp.tsv",
        diffexp_vs_dmrs_promoter="results/platforms_combined/varlo/rna_seq_comp/diffexp_vs_dmrs_promoter.tsv",
        gene_sets=config["enrichment"]["fgsea"]["gene_sets_file"],
        common_src=workflow.source_path("../scripts/common.R"),
    output:
        enrichment=report(
            "results/platforms_combined/varlo/rna_seq_comp/all-gene-sets.tsv",
            # caption="../report/fgsea-table-all.rst",
            category="DiffExp-Methylation Comparison",
            subcategory=f"Gene Set Enrichment Analysis",
            
        ),
        rank_ties=report(
            "results/platforms_combined/varlo/rna_seq_comp/rank-ties.tsv",
            # caption="../report/fgsea-rank-ties.rst",
            category="DiffExp-Methylation Comparison",
            subcategory=f"Gene Set Enrichment Analysis",
            
        ),
        significant=report(
            "results/platforms_combined/varlo/rna_seq_comp/sig-gene-sets.tsv",
            # caption="../report/fgsea-table-significant.rst",
            category="DiffExp-Methylation Comparison",
            subcategory=f"Gene Set Enrichment Analysis",
            
        ),
        plot=report(
            "results/platforms_combined/varlo/rna_seq_comp/table-plot.pdf",
            # caption="../report/fgsea-table-plot.rst",
            category="DiffExp-Methylation Comparison",
            subcategory=f"Gene Set Enrichment Analysis",
            
        ),
        plot_collapsed=report(
            "results/platforms_combined/varlo/rna_seq_comp/collapsed_pathways.table-plot.pdf",
            # caption="../report/fgsea-collapsed-table-plot.rst",
            category="DiffExp-Methylation Comparison",
            subcategory=f"Gene Set Enrichment Analysis",
            
        ),
    params:
        bioc_species_pkg=bioc_species_pkg,
        # model=get_model,
        gene_set_fdr=config["enrichment"]["fgsea"]["fdr_gene_set"],
        eps=config["enrichment"]["fgsea"]["eps"],
        # covariate=lambda w: config["diffexp"]["models"][w.model]["primary_variable"],
    conda:
        enrichment_env
    log:
        "logs/tables/fgsea/gene-set-enrichment.log",
    threads: 25
    script:
        "../scripts/fgsea.R"


# rule fgsea_plot_gene_sets:
#     input:
#         samples="results/sleuth/{model}.samples.tsv",
#         diffexp="results/tables/diffexp/{model}.genes-representative.diffexp.tsv",
#         gene_sets=config["enrichment"]["fgsea"]["gene_sets_file"],
#         sig_gene_sets="results/tables/fgsea/{model}.sig-gene-sets.tsv",
#         common_src=workflow.source_path("../scripts/common.R"),
#     output:
#         report(
#             directory("results/plots/fgsea/{model}"),
#             patterns=["{model}.{gene_set}.gene-set-plot.pdf"],
#             caption="../report/plot-fgsea-gene-set.rst",
#             category="Gene set enrichment analysis",
#             labels={"model": "{model}"},
#         ),
#     params:
#         model=get_model,
#         covariate=lambda w: config["diffexp"]["models"][w.model]["primary_variable"],
#     conda:
#         enrichment_env
#     log:
#         "logs/plots/fgsea/{model}.plot_fgsea_gene_set.log",
#     script:
#         "../scripts/plot-fgsea-gene-sets.R"


rule fgsea_datavzrd:
    input:
        config=workflow.source_path("../resources/fgsea.yaml"),
        comp="results/platforms_combined/varlo/rna_seq_comp/all-gene-sets.tsv",
    output:
        report(
            directory(
                "results/{platform}/{caller}/rna_seq_comp/gene_set"
            ),
            caption="../report/diffexp_vs_dmrs.rst",
            htmlindex="index.html",
            category="DiffExp-Methylation Comparison",
            subcategory="Gene Set Enrichment Analysis - table",
        ),
    log:
        "logs/diffexp_dmvzrd/diffexp_dmr_datavzrd/{platform}_{caller}.log",
    wrapper:
        "v3.13.0/utils/datavzrd"

