import sys

import pybedtools

sys.stderr = open(snakemake.log[0], "w", buffering=1)


def load_bedtool(path, name):
    """Load a BED file with better error handling."""
    try:
        return pybedtools.BedTool(path)
    except IndexError as e:
        print(f"Error loading {name} from {path}:", file=sys.stderr)
        print(
            f"  The file contains malformed BED records (missing required fields).",
            file=sys.stderr,
        )
        # Print the problematic lines
        try:
            with open(path, "r") as f:
                for i, line in enumerate(f, 1):
                    fields = line.strip().split("\t")
                    if len(fields) < 3:
                        print(
                            f"  Line {i}: {line.rstrip()} ({len(fields)} fields, expected ≥3)",
                            file=sys.stderr,
                        )
        except Exception as print_error:
            print(
                f"  Could not read file for diagnosis: {print_error}", file=sys.stderr
            )
        raise


try:
    this = load_bedtool(snakemake.input["this"], "this")
    other1 = load_bedtool(snakemake.input["other1"][0], "other1")
    other2 = load_bedtool(snakemake.input["other2"][0], "other2")

    # Find exclusive regions
    this_only = this.intersect(other1, v=True).intersect(other2, v=True)

    # Write output with validation
    valid_count = 0
    invalid_count = 0

    with open(snakemake.output[0], "w") as out:
        for region in this_only:
            mean_diff = float(region[4])
            # mean_methylation_difference should be between -1 and +1
            if mean_diff >= -1 and mean_diff <= 1:
                try:
                    print(region, file=out, end="")
                    valid_count += 1
                except (IndexError, ValueError):
                    invalid_count += 1
                    continue
            else:
                invalid_count += 1

    print(
        f"Wrote {valid_count} valid DMRs, filtered {invalid_count} invalid DMRs",
        file=sys.stderr,
    )

except Exception as e:
    print(f"Fatal error: {e}", file=sys.stderr)
    sys.exit(1)
