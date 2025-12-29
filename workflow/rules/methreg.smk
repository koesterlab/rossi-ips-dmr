# - Doku: https://www.bioconductor.org/packages/release/bioc/vignettes/MethReg/inst/doc/MethReg.html
# - GitHub: https://github.com/TransBioInfoLab/MethReg


rule methreg:
    input:
        # "results/methreg_data/clinical.rda",
        gene_exp_official="results/methreg_data/gene.exp.chr21.log2.rda",
        # All values in genes-aggregated are NA
        gene_exp="results/tables/diffexp/condition.genes-representative.diffexp.tsv",
        meth="results/{platform}/{caller}/meth_calling/calls.parquet",
        # calls="results/{platform}/{caller}/meth_calling/calls.parquet",
        meth_official="results/methreg_data/dna.met.chr21.rda",
    output:
        "results/{platform}/{caller}/methreg/methreg_results.parquet",
    conda:
        "../envs/methreg.yaml"
    log:
        "logs/methreg/{platform}_{caller}.log",
    script:
        "../scripts/methreg.R"
