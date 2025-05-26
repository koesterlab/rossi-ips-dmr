import pandas as pd

sys.stderr = open(snakemake.log[0], "w")

# Pfad zur Excel-Datei
rna_path = snakemake.input["rna_seq"]
genes_df = snakemake.input["genes_transcripts"]
germ_layer = snakemake.params["germ_layer"]

# Alle verfügbaren Tabellenblätter anzeigen
xls = pd.ExcelFile(rna_path)

# Ein bestimmtes Tabellenblatt einlesen (z. B. erstes Sheet)
rna_df = pd.read_excel(xls, sheet_name=xls.sheet_names[0])  # alternativ z. B. "Sheet1"
rna_df = genes_df[genes_df["test"] == germ_layer]
# Zeige die ersten 5 Zeilen zur Kontrolle
print(genes_df.head(), file=sys.stderr)
print(rna_df.head(), file=sys.stderr)
