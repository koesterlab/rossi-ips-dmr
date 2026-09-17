rule compute_own_biomarker:
    input:
        nanopore="results/nanopore/{caller}/meth_calling/calls_0.01.parquet",
        pacbio="results/pacbio/{caller}/meth_calling/calls_0.01.parquet",
    output:
        "results/platforms_combined/{caller}/plots_paper/own_biomarker.csv",
    resources:
        mem_mb=16000,
    conda:
        "../envs/python.yaml"
    log:
        "logs/compute_own_biomarker/{caller}.log",
    script:
        "../scripts/compute_own_biomarker.py"
