#### Used for debugging, not in use right now

import pandas as pd
import altair as alt

# Example for reading a BedGraph file without a header
file_path = snakemake.input[0]

# Manually adding the column names
column_names = [
    "chr",
    "start",
    "stop",
    "q-value",
    "mean_methylation_difference",
    "num_CpGs",
    "p_MWU",
    "p_2D_KS",
    "mean_g1",
    "mean_g2",
]
data = pd.read_csv(file_path, sep="\t", header=None, names=column_names)

# Check the data
print(data.head())

# Create the histogram for the column 'p_MWU'
hist_mwu = (
    alt.Chart(data)
    .mark_line()
    .encode(
        alt.X("p_MWU:Q", bin=alt.Bin(maxbins=30), title="p (MWU)"),
        alt.Y("count()", title="Count"),
    )
    .properties(title="Histogram of p (MWU)", width=300, height=200)
)

# Create the histogram for the column 'p_2D_KS'
hist_ks = (
    alt.Chart(data)
    .mark_line()
    .encode(
        alt.X("p_2D_KS:Q", bin=alt.Bin(maxbins=30), title="p (2D KS)"),
        alt.Y("count()", title="Count"),
    )
    .properties(title="Histogram of p (2D KS)", width=300, height=200)
)

# Display the histograms side by side
hist_mwu | hist_ks

# Combine the histograms in a vertical layout
combined_histograms = alt.vconcat(hist_mwu, hist_ks).resolve_scale(y="independent")

# Save the combined histograms as a PNG file
combined_histograms.save(snakemake.output[0], scale_factor=2.0)
