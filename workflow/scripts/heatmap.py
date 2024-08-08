import pandas as pd
import altair as alt
import os
from sklearn import cluster

filename_to_name = {
    "distal_intergenic": "Distal Intergenic",
    "promoter": "Promoter",
    "intron": "Intron",
    "exon": "Exon",
    "3_utr": "3' UTR",
    "5_utr": "5' UTR",
    "downstream": "Downstream"
}

# Group df by gene regions
def aggregate_by_gene_region(df):
    df['transcriptId'] = df['transcriptId'].astype(str)
    df['annotation'] = df['annotation'].astype(str)
    # For every gene region take the maximum of methylation differences
    df_grouped = df.groupby(['transcriptId', 'annotation'])[
        'mean_methylation_difference'].max().reset_index()
    df_grouped['annotation_type'] = df_grouped['annotation'].str.replace("\s*\([^)]*\)", "", regex=True)
    df_grouped['region'] = df_grouped.apply(
        lambda row: f"{row['transcriptId']}:{row['annotation']}", axis=1)    
    return df_grouped[['region', 'annotation_type', 'mean_methylation_difference']]

def cluster_data(df):
    # Perform clustering
    df_cluster = df.drop(['region', 'annotation_type'], axis=1)
    # 5 Cluster because: All same, 1&2, 2&3, 1&3, all different
    # Or 2 Clusters because SSE is best there
    k_means = cluster.KMeans(n_clusters=2)
    k_means.fit(df_cluster) 
    labels = k_means.labels_
    heatmap_data['cluster'] = labels
    return df.sort_values('cluster')
    # print(heatmap_data.to_string())

input_files = snakemake.input
output = snakemake.output[0]

# Merge single dfs to one df
aggregated_data = [aggregate_by_gene_region(
    pd.read_csv(file, sep='\t')) for file in input_files]
all_regions = set(region for df in aggregated_data for region in df['region'])
for idx, df in enumerate(aggregated_data):
    aggregated_data[idx] = df.rename(columns={"mean_methylation_difference": f'BC0{idx+2}'})
heatmap_data = aggregated_data[0].merge(aggregated_data[1], on=['region', 'annotation_type'], how='outer').merge(aggregated_data[2], on=['region', 'annotation_type'], how='outer')
heatmap_data = heatmap_data.fillna(0)


heatmap_data = cluster_data(heatmap_data)

# import matplotlib.pyplot as plt
# numClusters = [1,2,3,4,5,6, 7, 8, 9, 10]
# SSE = []
# for k in numClusters:
#     k_means = cluster.KMeans(n_clusters=k)
#     k_means.fit(df_cluster)
#     SSE.append(k_means.inertia_)

# plt.plot(numClusters, SSE)
# plt.xlabel('Number of Clusters')
# plt.ylabel('SSE')
# plt.show()


# Convert to plotable form
df_melted = heatmap_data.melt(id_vars=['region', 'annotation_type', 'cluster'], 
                    value_vars=['BC02', 'BC03', 'BC04'], 
                    var_name='sample', 
                    value_name='Methylation difference')

# Filter on annotation_type
name = os.path.basename(output).replace('.png', '')
annotation_type = filename_to_name[name]
df_filtered = df_melted[df_melted['annotation_type'] == annotation_type]

# Create heatmap
color_scale = alt.Scale(domain=[-1, 0, 1], range=['red', 'white', 'blue'])
heatmap = alt.Chart(df_filtered).mark_rect().encode(
    x=alt.X('sample:N', title='Samples'),
    y=alt.Y('region:N', title='Generegions', axis=alt.Axis(labels=False)),
    color=alt.Color('Methylation difference:Q',
                    scale=color_scale, title='Methylation Difference')
).properties(
    title=f'DMRs in {annotation_type}',
    width=800,
    height=600
)

heatmap.save(output, format='png')
