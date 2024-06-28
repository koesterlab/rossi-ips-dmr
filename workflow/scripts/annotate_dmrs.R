# Pakete installieren, falls noch nicht geschehen
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("genomation", "methylKit"))

# Bibliotheken laden
library(genomation)
library(methylKit)

# DMR-Daten aus Datei einlesen
dmr_file <- "/projects/koesterlab/benchmark-methylation/data_analysis_jochen/compare_meth/results/dmr_calls/metilene_output_BC02-ref-sorted.txt"
dmr_data <- read.table(dmr_file, header = TRUE)

# DMR-Daten in ein GRanges-Objekt konvertieren
dmr_gr <- GRanges(
  seqnames = Rle(dmr_data$chr),
  ranges = IRanges(start = dmr_data$start, end = dmr_data$stop),
  q_value = dmr_data$q.value,
  mean_meth_diff = dmr_data$mean.methylation.difference,
  CpGs = dmr_data$X.CpGs,
  p_MWU = dmr_data$p..MWU.,
  p_2D_KS = dmr_data$p..2D.KS.,
  mean_g1 = dmr_data$mean.g1,
  mean_g2 = dmr_data$mean.g2
)

# Genom-Annotationsdatei (z.B. BED-Datei) einlesen
transcriptBED <- system.file("extdata", "/projects/koesterlab/benchmark-methylation/data_analysis_jochen/compare_meth/results/dmr_calls/metilene_output_BC02-ref-sorted.txt", package = "genomation")
gene_obj <- readTranscriptFeatures(transcriptBED)

# DMRs annotieren mit Promotor/Exon/Intron-Informationen
annotated_dmrs <- annotateWithGeneParts(dmr_gr, gene_obj)

# Ergebnisse anzeigen
print(annotated_dmrs)

# CpG Inseln und Umlandregionen (Shores) einlesen und annotieren
cpg_file <- system.file("extdata", "/projects/koesterlab/benchmark-methylation/data_analysis_jochen/compare_meth/results/dmr_calls/metilene_output_BC02-ref-sorted.txt", package = "genomation")
cpg_obj <- readFeatureFlank(cpg_file, feature.flank.name = c("CpGi", "shores"))

# DMRs annotieren mit CpG Inseln und Shores
diffCpGann <- annotateWithFeatureFlank(
  dmr_gr,
  cpg_obj$CpGi,
  cpg_obj$shores,
  feature.name = "CpGi",
  flank.name = "shores"
)

# Ergebnisse anzeigen
print(diffCpGann)
