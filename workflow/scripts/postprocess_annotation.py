import pandas as pd
import numpy as np

# Read the TSV file
df = pd.read_csv(snakemake.input[0], sep='\t')

# Function to split the attributes column into separate columns


def split_attributes(attributes):
    # Split the attributes string by semicolon
    attr_dict = {}
    for item in attributes.split(';'):
        if '=' in item:
            key, value = item.split('=')
            attr_dict[key] = value
    return attr_dict


# Apply the function to the attributes column
attributes_df = df['attributes'].apply(split_attributes).apply(pd.Series)

# Merge the new columns with the original dataframe
df = pd.concat([df.drop(columns=['attributes']), attributes_df], axis=1)

# Remove rows where 'Name' column is NaN
if 'Name' in df.columns:
    df = df.dropna(subset=['Name'])

if 'gene_name' in df.columns:
    df = df.dropna(subset=["gene_name"])

# Vertausche die Reihenfolge der Spalten gene_name und gene_biotype
if 'gene_name' in df.columns and 'gene_biotype' in df.columns:
    columns = df.columns.tolist()
    gene_name_index = columns.index('gene_name')
    gene_biotype_index = columns.index('gene_biotype')
    # Vertausche die Positionen der beiden Spalten
    columns[gene_name_index], columns[gene_biotype_index] = columns[gene_biotype_index], columns[gene_name_index]
    df = df[columns]

# Remove rows where 'type' column is either 'chromosome' or 'exon'
df = df[~df['type'].isin(['chromosome', 'exon'])]


# Remove rows where 'ID' column starts with 'gene'
df = df[~df['ID'].fillna('').astype(str).str.startswith('gene')]

if 'Parent' in df.columns:
    df['Parent'] = df['Parent'].str.replace('gene:', '', regex=True)

df = df[df['q-value'] <= 0.25]

# Sort dataframe by 'q-value' in ascending order

# Remove rows where 'q-value' > 0.25
# Berechnung der neuen Spalte absolute_signed_pi_val
df['absolute_signed_pi_val'] = abs(-np.log10(df['p(MWU)'])
                                   * df['mean_methylation_difference'])

# Sortierung der Tabelle aufsteigend nach der neuen Spalte
df = df.sort_values(by='absolute_signed_pi_val', ascending=True)


df = df.dropna(axis=1, how='all')
# Save the new dataframe to a TSV file
df.to_csv(snakemake.output[0], sep='\t', index=False)
