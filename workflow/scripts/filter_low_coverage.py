sys.stderr = open(snakemake.log[0], "w", buffering=1)

input_vcf = snakemake.input[0]  
output_vcf = snakemake.output[0]

with open(input_vcf, "r") as infile, open(output_vcf, "w") as outfile:
    for line in infile:
        if line.startswith("##"):
            outfile.write(line)
        elif line.startswith("#"):
            outfile.write(line)
        else:
            fields = line.strip().split("\t")
            info_col = fields[8].split(":")
            sample_col = fields[9].split(":")
            
            if "DP" in info_col:
                dp_index = info_col.index("DP")
                dp_value = int(sample_col[dp_index])
                
                if dp_value > 20:
                    outfile.write(line)