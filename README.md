# Snakemake workflow: rossi-ips-dmr

A Snakemake workflow that compares DNA methylation and gene expression of human
induced pluripotent stem cells (iPSCs, `psc`) and their differentiated germ layers
(`endoderm`, `mesoderm`, `ectoderm`).

## Overview

1. **Methylation calling** with [varlociraptor](https://varlociraptor.github.io) on
   PacBio and Nanopore alignments, per platform and jointly (`platforms_combined`).
   Candidates are split into `scatter_items` chunks and gathered again.
2. **DMR calling** with [metilene](http://www.bioinf.uni-leipzig.de/Software/metilene/)
   between a base experiment (`ACTIVE_BASE_COMP` in `workflow/rules/common.smk`,
   currently `psc`) and each other germ layer. DMRs overlapping DMRs of the other
   layers are removed (layer-specific DMRs).
3. **Annotation** of DMRs with genomic elements (ChIPseeker), gene names (from the
   Ensembl GTF) and Ensembl regulatory features.
4. **Differential expression** of long-read RNA-seq data via the
   [rna-seq-kallisto-sleuth](https://github.com/snakemake-workflows/rna-seq-kallisto-sleuth)
   workflow, used as a Snakemake module.
5. **Comparison** of methylation differences and differential expression, plus
   plots for the paper, all collected in a Snakemake report.

## Input

- Alignments with methylation tags: `resources/{platform}/{germ_layer}.bam`
  (`platform` is `pacbio` or `nanopore`).
- RNA-seq data: SRA accessions (`rna_accessions_old`) and a zip archive with
  Nanopore BAMs at `resources/rna_seq_new/KOLF_Trilineage_RNAseq_new.zip`
  (`rna_accessions_new`). `config/samples.tsv` and `config/units.tsv` for the
  kallisto-sleuth module are generated automatically.
- Reference genome, annotation and regulatory features are downloaded from Ensembl
  (release in `resources: ref:` of `config/config.yaml`).

## Configuration

All settings are in `config/config.yaml`.

## Usage

```bash
snakemake --sdm conda --cores <n>
snakemake --report report.zip
```
The usage of this workflow is described in the [Snakemake Workflow Catalog](https://snakemake.github.io/snakemake-workflow-catalog/?usage=<owner>%2F<repo>).

If you use this workflow in a paper, don't forget to give credits to the authors by citing the URL of this (original) <repo>sitory and its DOI (see above).
