rule rna_seq:
    input:
        genes_transcripts="results/{platform}/{caller}/dmr_calls/{group2}/genes_transcripts/chipseeker_postprocessed.tsv",
        rna_seq="resources/rna_seq/rna_seq.xlsx",
    output:
        "results/{platform}/{caller}/rna_seq/{group2}.csv",
    params:
        groupt2=lambda wildcards: wildcards.group2,
    log:
        "logs/rna_seq/{platform}/{caller}/{group2}.log",
    script:
        "../scripts/compare_rna_seq.py"
