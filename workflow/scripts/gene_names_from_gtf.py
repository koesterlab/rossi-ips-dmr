import sys

import polars as pl

sys.stderr = open(snakemake.log[0], "w", buffering=1)

transcript_ids = (
    pl.read_csv(snakemake.input.chipseeker, separator="\t", columns=["transcriptId"])
    .get_column("transcriptId")
    .unique()
)

# GTF columns: seqname, source, feature, start, end, score, strand, frame, attributes
gene_names = (
    pl.read_csv(
        snakemake.input.gtf,
        separator="\t",
        has_header=False,
        comment_prefix="#",
        quote_char=None,
        columns=[2, 8],
        new_columns=["feature", "attributes"],
    )
    .filter(pl.col("feature") == "transcript")
    .select(
        pl.col("attributes")
        .str.extract(r'transcript_id "([^"]+)"')
        .alias("ensembl_transcript_id"),
        pl.col("attributes").str.extract(r'gene_id "([^"]+)"').alias("ensembl_gene_id"),
        # Not every gene has a name; missing names are written as empty strings
        pl.col("attributes")
        .str.extract(r'gene_name "([^"]+)"')
        .alias("external_gene_name"),
    )
    .filter(pl.col("ensembl_transcript_id").is_in(transcript_ids.implode()))
)

print(
    f"Found {gene_names.height} of {len(transcript_ids)} transcripts in the annotation",
    file=sys.stderr,
)
# Same columns as the previous biomaRt query, so downstream scripts stay unchanged
gene_names.write_csv(snakemake.output[0], separator="\t")
