import sys

import decoupler as dc

sys.stderr = open(snakemake.log[0], "w", buffering=1)

# TF-target interactions (columns source, target, weight) from CollecTRI
net = dc.op.collectri(organism="human")
net.to_csv(snakemake.output[0])
