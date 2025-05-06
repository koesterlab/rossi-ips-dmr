import pysam

# Pfad zur BCF-Datei
bcf_file = snakemake.input[0]

# BCF-Datei öffnen
bcf = pysam.VariantFile(bcf_file)  # pysam erkennt das Format automatisch

# BED-Datei erstellen
with open(snakemake.output[0], "w") as bed_file:
    for record in bcf:
        # Extrahiere benötigte Informationen aus der BCF-Datei
        chrom = record.chrom
        start = record.pos - 1  # BED ist 0-basiert
        end = record.pos  # VCF/BCF ist 1-basiert

        # Schreibe die Daten im BED-Format
        bed_file.write(f"{chrom}\t{start}\t{end}\n")
