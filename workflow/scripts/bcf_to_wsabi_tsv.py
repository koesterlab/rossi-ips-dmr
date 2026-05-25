import gzip
import sys

import pysam

sys.stderr = open(snakemake.log[0], "w", buffering=1)

vcf = pysam.VariantFile(snakemake.input[0])

with gzip.open(snakemake.output[0], "wt") as outf:
    outf.write("C\tP\tB\n")

    for record in vcf:
        if record.alts and record.alts[0] == "<METH>":
            chrom = record.chrom
            pos = record.pos

            af_values = []
            for sample in record.samples.values():
                try:
                    af_values.append(float(sample.get("AF", [0])[0]))
                except (ValueError, TypeError):
                    pass

            # Average the AF values across samples
            if af_values:
                avg_af = sum(af_values) / len(af_values)
                outf.write(f"{chrom}\t{pos}\t{avg_af}\n")

vcf.close()
