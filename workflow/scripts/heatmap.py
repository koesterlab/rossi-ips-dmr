import pandas as pd
import altair as alt

# Group df by gene regions
def aggregate_by_gene_region(df):
    df['geneChr'] = df['geneChr'].astype(str)
    df['geneStart'] = df['geneStart'].astype(str)
    df['geneEnd'] = df['geneEnd'].astype(str)
    # For every gene region compute mean of methylation difference
    df_grouped = df.groupby(['geneChr', 'geneStart', 'geneEnd'])['mean_methylation_difference'].mean().reset_index()
    df_grouped['region'] = df_grouped.apply(lambda row: f"{row['geneChr']}:{row['geneStart']}-{row['geneEnd']}", axis=1)
    return df_grouped[['region', 'mean_methylation_difference']]


input_files = snakemake.input
output_file = snakemake.output[0]

aggregated_data = [aggregate_by_gene_region(pd.read_csv(file, sep='\t')) for file in input_files]

all_regions = set(region for df in aggregated_data for region in df['region'])
heatmap_data = pd.DataFrame(index=sorted(all_regions))

# Create one dataset containing all samples
for idx, df in enumerate(aggregated_data):
    # Make regions to index of df
    df_pivot = df.pivot_table(index='region', values='mean_methylation_difference')
    heatmap_data[f'BC0{idx+2}'] = heatmap_data.index.map(df_pivot['mean_methylation_difference'].to_dict())

heatmap_data = heatmap_data.fillna(0).reset_index().melt(id_vars='index', var_name='Sample', value_name='Mean_Methylation_Difference')

# Create heatmap
color_scale = alt.Scale(domain=[-1, 0, 1], range=['red', 'white', 'blue'])

heatmap = alt.Chart(heatmap_data).mark_rect().encode(
    x=alt.X('Sample:N', title='Samples'),
    y=alt.Y('index:N', title='Genregionen'),
    color=alt.Color('Mean_Methylation_Difference:Q', scale=color_scale, title='Mean Methylation Difference')
).properties(
    title='DMRs relative to next gene',
    width=800,
    height=600
)

heatmap.save(output_file, format='png')